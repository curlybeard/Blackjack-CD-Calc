#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <unordered_map>
#include <tuple>
#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <chrono>
#include <cstring>
#include <limits>
#include <cctype>
#include <cstdarg>
#include <sstream>
#include <omp.h> 

using namespace std;

/*
COMPILE:
  g++ -O3 -fopenmp blackjackCDC.cpp -o blackjackCDC
*/


//define all possible card 11 integer face card values in a deck (1-10)
typedef array<int, 11> Deck;

//hash deck to find identical deck compositions to skip computing duplicate deck compositions
struct DeckHash 
{
    size_t operator()(const Deck& d) const 
    {
        size_t h = 0;
        for (int i = 0; i < 11; i++) 
        {
            h ^= std::hash<int>{}(d[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};


unordered_map<Deck, int, DeckHash> indexer;
unordered_map<Deck, int, DeckHash> cache_idx;

//number of unique deck states
int deck_number = 1;
//amount of memory allocated for theoretical max number of deck and dealer compositions with upper bound of 8 deck and 50% pen
const int MAX_DECKS = 40000; 
const int MAX_COMPS = 3200;

//pointer-based flat array faster that a map
int* deck_mapping_flat;
int* decks_flat;

double* cache_dealer_flat;
double* cache_stand_flat;
double* cache_hit_hard_flat;
double* cache_hit_soft_flat;
double* cache_double_hard_flat;
double* cache_double_soft_flat;
double* cache_split_flat;

int composition_counts[11][MAX_COMPS];
int composition_values[11][MAX_COMPS];
int composition_cards[MAX_COMPS][11];
vector<int> composition_lists[11][6];
int comp_number = 1;

//in case max_decks exceeded, but shouldn't be
static inline void fatalf(const char* fmt, ...) 
{
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
    exit(1);
}

//ensure card face values are strictly 1-10
static inline void require_card_val(int c, const char* label) 
{
    if (c < 1 || c > 10) fatalf("Invalid %s: %d (must be 1..10)", label, c);
}

static inline double qnan() { return std::numeric_limits<double>::quiet_NaN(); }

static inline void fill_nan(double* p, size_t n) 
{
    const double nanv = qnan();
    for (size_t i = 0; i < n; i++) p[i] = nanv;
}

//hybrid ENHC, dealer only checks for BJ on Ace, not 10
static inline bool dealer_peeks_for_bj(int upcard) 
{
    return upcard == 1;
}

//strucut to hold all possible player decisions
struct ActionEVs 
{
    double stand = NAN, hit = NAN, dbl = NAN, split = NAN, dealer_bj_prob_now = 0.0;
    bool player_blackjack = false;
};


//compare EV on player actions, and return highest EV action
static inline const char* best_action_name(const ActionEVs& a, double& best_ev_out) 
{
    const char* best_name = "Stand";
    double best = a.stand;

    auto consider = [&](const char* nm, double v) 
    {
        if (!std::isnan(v) && (std::isnan(best) || v > best)) 
        {
            best = v;
            best_name = nm;
        }
    };

    consider("Hit", a.hit);
    consider("Double", a.dbl);
    consider("Split", a.split);

    best_ev_out = best;
    return best_name;
}

//return numerical EV value of highest EV action
static inline double best_ev_value(const ActionEVs& a)
{
    bool any = false;
    double best = -1e300;

    auto take = [&](double v) 
    {
        if (std::isfinite(v)) 
        {
            if (!any || v > best) best = v;
            any = true;
        }
    };

    take(a.stand);
    take(a.hit);
    take(a.dbl);
    take(a.split);

    return any ? best : NAN;
}

//pre-allocate memroy for caches to prevent OS memory allocation mid calculation
void initialize() 
{
    deck_mapping_flat = new int[MAX_DECKS * 11]();
    decks_flat        = new int[MAX_DECKS * 11]();

    cache_dealer_flat      = new double[(size_t)MAX_DECKS * 10 * 6];
    cache_stand_flat       = new double[(size_t)MAX_DECKS * 10 * 6];
    cache_hit_hard_flat    = new double[(size_t)MAX_DECKS * 10 * 18];
    cache_hit_soft_flat    = new double[(size_t)MAX_DECKS * 10 * 10];
    cache_double_hard_flat = new double[(size_t)MAX_DECKS * 10 * 18];
    cache_double_soft_flat = new double[(size_t)MAX_DECKS * 10 * 10];
    cache_split_flat       = new double[(size_t)MAX_DECKS * 10 * 10];

    fill_nan(cache_dealer_flat,      (size_t)MAX_DECKS * 10 * 6);
    fill_nan(cache_stand_flat,       (size_t)MAX_DECKS * 10 * 6);
    fill_nan(cache_hit_hard_flat,    (size_t)MAX_DECKS * 10 * 18);
    fill_nan(cache_hit_soft_flat,    (size_t)MAX_DECKS * 10 * 10);
    fill_nan(cache_double_hard_flat, (size_t)MAX_DECKS * 10 * 18);
    fill_nan(cache_double_soft_flat, (size_t)MAX_DECKS * 10 * 10);
    fill_nan(cache_split_flat,       (size_t)MAX_DECKS * 10 * 10);

    memset(composition_counts, 0, sizeof(composition_counts));
    memset(composition_values, 0, sizeof(composition_values));
    memset(composition_cards,  0, sizeof(composition_cards));
}

//simulate dealer placing a card into player/dealer hands
inline Deck add(Deck d, int c) { d[0]++; d[c]++; return d; }
//remove card form remaining shoe
inline Deck draw(Deck d, int c) { d[0]--; d[c]--; return d; }

//recursively calculate outcome of all dealer hands (17,18,19,20,21,bust)
int compute_dealer_sequences(int first, Deck cards, int total, bool soft = false) 
{
    if (total > 21 && soft) return compute_dealer_sequences(first, cards, total - 10, false);
    if (total >= 17) 
    {
        static unordered_map<Deck, int, DeckHash> comp_indexer;
        if (comp_indexer[cards] == 0) 
        {
            comp_indexer[cards] = comp_number++;
            if (comp_number >= MAX_COMPS) fatalf("MAX_COMPS exceeded");
            for (int i = 0; i < 11; i++) composition_cards[comp_indexer[cards]][i] = cards[i];
        }
        int idx = comp_indexer[cards];
        composition_values[first][idx] = min(total, 22);
        composition_counts[first][idx]++;
        return 1;
    }
    int ans = 0;
    for (int c = 1; c < 11; c++) 
    {
        Deck c2 = add(cards, c);
        if (soft && c == 1) ans += compute_dealer_sequences(first, c2, total + 1, true);
        else if (c == 1)    ans += compute_dealer_sequences(first, c2, total + 11, true);
        else                ans += compute_dealer_sequences(first, c2, total + c, soft);
    }
    return ans;
}

//sort dealer outcomes based on dealer upcard and total value
void dealer_sequences() 
{
    Deck empty_deck = {0};
    for (int c = 1; c < 11; c++) compute_dealer_sequences(c, empty_deck, (c == 1 ? 11 : c), c == 1);

    for (int c = 1; c < 11; c++) 
    {
        for (int comp = 1; comp < comp_number; comp++) 
        {
            if (composition_counts[c][comp] > 0) 
            {
                composition_lists[c][composition_values[c][comp] - 17].push_back(comp);
            }
        }
    }
}

//calculate probability of dealer getting a natural blackjack ie) Ace and 10 or 10 and Ace
inline double dealer_natural_probability(int deck_idx, int upcard) 
{
    int* d = &decks_flat[deck_idx * 11];
    if (upcard == 1)  return (d[0] > 0) ? ((double)d[10] / d[0]) : 0.0;
    if (upcard == 10) return (d[0] > 0) ? ((double)d[1] / d[0])  : 0.0;
    return 0.0;
}


//calculate exact chance of a dealer reaching 17,18,19,20,21 or busting given they do not have a natural blackjack
double dealer_conditional_probability(int deck, int upcard, int value) 
{
    int cache_off = (deck * 60) + ((upcard - 1) * 6) + (value - 17);
    if (!std::isnan(cache_dealer_flat[cache_off])) return cache_dealer_flat[cache_off];

    int* d = &decks_flat[deck * 11];
    double bj = dealer_natural_probability(deck, upcard);
    double ans = (value == 21) ? -bj : 0.0;

    for (int comp : composition_lists[upcard][value - 17]) 
    {
        int* cards = composition_cards[comp];
        double proba = composition_counts[upcard][comp];
        double denom = d[0];
        bool valid = true;
        for (int i = 1; i < 11; i++) 
        {
            if (d[i] >= cards[i]) 
            {
                for (int j = 0; j < cards[i]; j++) 
                {
                    proba *= (d[i] - j) / denom;
                    denom -= 1.0;
                }
            } 
            else 
            { 
                valid = false; break; 
            }
        }
        if (valid) ans += proba;
    }
    //divide by probability dealer doesn't have bj
    ans /= (1.0 - bj);
    return cache_dealer_flat[cache_off] = ans;
}

//calculate EV of the stand action 
double ev_stand(int deck, int total, int upcard) 
{
    int cache_val = max(total - 16, 0);
    int cache_off = (deck * 60) + ((upcard - 1) * 6) + cache_val;
    if (!std::isnan(cache_stand_flat[cache_off])) return cache_stand_flat[cache_off];

    const double bj = dealer_natural_probability(deck, upcard);
    const bool peek = dealer_peeks_for_bj(upcard);

    //if dealer peeks (Ace), decisions happen only after no BJ is confirmed
    //therefore, action EVs should be conditional on no dealer BJ, with no standalone bj term
    if (peek) 
    {
        if (bj >= 1.0 - 1e-15) return cache_stand_flat[cache_off] = -1.0;

        const double p17   = dealer_conditional_probability(deck, upcard, 17);
        const double p18   = dealer_conditional_probability(deck, upcard, 18);
        const double p19   = dealer_conditional_probability(deck, upcard, 19);
        const double p20   = dealer_conditional_probability(deck, upcard, 20);
        const double p21nb = dealer_conditional_probability(deck, upcard, 21);
        const double pbust = dealer_conditional_probability(deck, upcard, 22);

        double ev = pbust;
        if (17 < total) ev += p17; else if (17 > total) ev -= p17;
        if (18 < total) ev += p18; else if (18 > total) ev -= p18;
        if (19 < total) ev += p19; else if (19 > total) ev -= p19;
        if (20 < total) ev += p20; else if (20 > total) ev -= p20;
        if (21 < total) ev += p21nb; else if (21 > total) ev -= p21nb;

        return cache_stand_flat[cache_off] = ev;
    }

    //no peek (10 upcard hybrid case): ENHC behavior (unconditional EV with -bj term)
    if (bj >= 1.0 - 1e-15) return cache_stand_flat[cache_off] = -1.0;

    double ombj = 1.0 - bj;
    double p17   = ombj * dealer_conditional_probability(deck, upcard, 17);
    double p18   = ombj * dealer_conditional_probability(deck, upcard, 18);
    double p19   = ombj * dealer_conditional_probability(deck, upcard, 19);
    double p20   = ombj * dealer_conditional_probability(deck, upcard, 20);
    double p21nb = ombj * dealer_conditional_probability(deck, upcard, 21);
    double pbust = ombj * dealer_conditional_probability(deck, upcard, 22);

    double ev = -bj + pbust;
    if (17 < total) ev += p17; else if (17 > total) ev -= p17;
    if (18 < total) ev += p18; else if (18 > total) ev -= p18;
    if (19 < total) ev += p19; else if (19 > total) ev -= p19;
    if (20 < total) ev += p20; else if (20 > total) ev -= p20;
    if (21 < total) ev += p21nb; else if (21 > total) ev -= p21nb;

    return cache_stand_flat[cache_off] = ev;
}

double ev_hit_soft(int deck, int total, int upcard);
double ev_hit_hard(int deck, int total, int upcard);

//resurisvely calculate EV of hit action for a hard total
double ev_hit_hard(int deck, int total, int upcard) 
{
    int cache_off = (deck * 180) + ((upcard - 1) * 18) + (total - 4);
    if (!std::isnan(cache_hit_hard_flat[cache_off])) return cache_hit_hard_flat[cache_off];

    int* d = &decks_flat[deck * 11];
    double ev = 0.0;
    for (int c = 1; c < 11; c++) 
    {
        if (d[c] > 0) 
        {
            int deck2 = deck_mapping_flat[deck * 11 + c];
            double p = (double)d[c] / d[0];
            if (c == 1 && total + 11 <= 21)
                ev += p * max(ev_hit_soft(deck2, total + 11, upcard), ev_stand(deck2, total + 11, upcard));
            else if (total + c <= 21)
                ev += p * max(ev_hit_hard(deck2, total + c, upcard), ev_stand(deck2, total + c, upcard));
            else
                ev -= p;
        }
    }
    return cache_hit_hard_flat[cache_off] = ev;
}

//resursively calculate EV of hit action for a soft total
double ev_hit_soft(int deck, int total, int upcard) 
{
    int cache_off = (deck * 100) + ((upcard - 1) * 10) + (total - 12);
    if (!std::isnan(cache_hit_soft_flat[cache_off])) return cache_hit_soft_flat[cache_off];

    int* d = &decks_flat[deck * 11];
    double ev = 0.0;
    for (int c = 1; c < 11; c++) 
    {
        if (d[c] > 0) 
        {
            double p = (double)d[c] / d[0];
            int deck2 = deck_mapping_flat[deck * 11 + c];
            if (total + c <= 21)
                ev += p * max(ev_hit_soft(deck2, total + c, upcard), ev_stand(deck2, total + c, upcard));
            else
                ev += p * max(ev_hit_hard(deck2, total + c - 10, upcard), ev_stand(deck2, total + c - 10, upcard));
        }
    }
    return cache_hit_soft_flat[cache_off] = ev;
}

//calculate EV of double action of a hard total
double ev_double_hard(int deck, int total, int upcard) 
{
    int cache_off = (deck * 180) + ((upcard - 1) * 18) + (total - 4);
    if (!std::isnan(cache_double_hard_flat[cache_off])) return cache_double_hard_flat[cache_off];

    int* d = &decks_flat[deck * 11];
    double ev = 0.0;
    for (int c = 1; c < 11; c++) 
    {
        if (d[c] > 0) 
        {
            int deck2 = deck_mapping_flat[deck * 11 + c];
            double p = (double)d[c] / d[0];
            if (c == 1 && total + 11 <= 21) ev += p * ev_stand(deck2, 11 + total, upcard);
            else if (total + c <= 21)       ev += p * ev_stand(deck2, total + c, upcard);
            else                            ev -= p;
        }
    }
    return cache_double_hard_flat[cache_off] = 2 * ev;
}

// calculate EV of double action of a soft total
double ev_double_soft(int deck, int total, int upcard) 
{
    int cache_off = (deck * 100) + ((upcard - 1) * 10) + (total - 12);
    if (!std::isnan(cache_double_soft_flat[cache_off])) return cache_double_soft_flat[cache_off];

    int* d = &decks_flat[deck * 11];
    double ev = 0.0;
    for (int c = 1; c < 11; c++) {
        if (d[c] > 0) {
            double p = (double)d[c] / d[0];
            int deck2 = deck_mapping_flat[deck * 11 + c];
            if (total + c <= 21) ev += p * ev_stand(deck2, total + c, upcard);
            else                 ev += p * ev_stand(deck2, total + c - 10, upcard);
        }
    }
    return cache_double_soft_flat[cache_off] = 2 * ev;
}

//calculate EV of split action, one split only, multiply ev of first hand by 2
double ev_split(int deck, int card, int upcard) 
{
    int cache_off = (deck * 100) + ((upcard - 1) * 10) + (card - 1);
    if (!std::isnan(cache_split_flat[cache_off])) return cache_split_flat[cache_off];

    int* d = &decks_flat[deck * 11];
    double ev = 0.0;
    for (int c = 1; c < 11; c++) 
    {
        if (d[c] > 0) 
        {
            double p = (double)d[c] / d[0];
            int deck2 = deck_mapping_flat[deck * 11 + c];
            //for split aces, only one card allowed
            if (card == 1) ev += p * ev_stand(deck2, 11 + c, upcard);
            else if (c == 1) 
            {
                int st = 11 + card;
                ev += p * max({ev_stand(deck2, st, upcard),
                               ev_hit_soft(deck2, st, upcard),
                               ev_double_soft(deck2, st, upcard)});
            } 
            else 
            {
                int ht = card + c;
                ev += p * max({ev_stand(deck2, ht, upcard),
                               ev_hit_hard(deck2, ht, upcard),
                               ev_double_hard(deck2, ht, upcard)});
            }
        }
    }
    return cache_split_flat[cache_off] = 2 * ev;
}

//store all possible EVs for a specific hand into ActionEV struct
static ActionEVs compute_action_evs(int deck_idx, int p1, int p2, int upcard) {
    ActionEVs out;

    const double bj_prob = dealer_natural_probability(deck_idx, upcard);
    const bool peek = dealer_peeks_for_bj(upcard);

    //after a peek on Ace, decisions only occur given "no dealer BJ"
    out.dealer_bj_prob_now = peek ? 0.0 : bj_prob;

    //player blackjack
    if ((p1 == 1 && p2 == 10) || (p1 == 10 && p2 == 1)) 
    {
        out.player_blackjack = true;
        out.stand = (1.0 - bj_prob) * 1.5;
        return out;
    }

    //compute conditional (peek) or unconditional (no-peek) action EVs via ev_* tree.
    if (p1 == p2) 
    {
        if (p1 == 1) 
        {
            out.hit   = ev_hit_soft(deck_idx, 12, upcard);
            out.stand = ev_stand(deck_idx, 12, upcard);
            out.dbl   = ev_double_soft(deck_idx, 12, upcard);
            out.split = ev_split(deck_idx, 1, upcard);
        } 
        else 
        {
            out.hit   = ev_hit_hard(deck_idx, p1 * 2, upcard);
            out.stand = ev_stand(deck_idx, p1 * 2, upcard);
            out.dbl   = ev_double_hard(deck_idx, p1 * 2, upcard);
            out.split = ev_split(deck_idx, p1, upcard);
        }
    } 
    else if (p1 == 1 || p2 == 1) 
    {
        int st = p1 + p2 + 10;
        out.hit   = ev_hit_soft(deck_idx, st, upcard);
        out.stand = ev_stand(deck_idx, st, upcard);
        out.dbl   = ev_double_soft(deck_idx, st, upcard);
    } 
    else 
    {
        int ht = p1 + p2;
        out.hit   = ev_hit_hard(deck_idx, ht, upcard);
        out.stand = ev_stand(deck_idx, ht, upcard);
        out.dbl   = ev_double_hard(deck_idx, ht, upcard);
    }

    // Hybrid adjustment for peek-on-Ace:
    // True EV from the start of the round is:
    //   EV = P(BJ)*(-1) + P(noBJ)*EV_conditional
    //
    // Crucially: the -1 loss does NOT scale with double/split because extra bets are never placed
    // if the dealer has BJ (peek ends the round immediately).
    if (peek) 
    {
        auto wrap = [&](double v) -> double 
        {
            if (!std::isfinite(v)) return v;
            return (-bj_prob) + (1.0 - bj_prob) * v;
        };
        out.stand = wrap(out.stand);
        out.hit   = wrap(out.hit);
        out.dbl   = wrap(out.dbl);
        out.split = wrap(out.split);
    }

    return out;
}

//recursively map all possible deck states to a unique index
int recurse_decks(Deck d, int total) 
{
    if (indexer.count(d) && cache_idx[d] <= total) return indexer[d];

    if (indexer[d] == 0) 
    {
        if (deck_number >= MAX_DECKS) fatalf("Fatal: MAX_DECKS (%d) exceeded.", MAX_DECKS);
        indexer[d] = deck_number++;
        for (int i = 0; i < 11; i++) decks_flat[indexer[d] * 11 + i] = d[i];
    }
    int id = indexer[d];
    cache_idx[d] = total;

    for (int c = 1; c < min(22 - total, 11); c++) 
    {
        if (d[c] > 0) deck_mapping_flat[id * 11 + c] = recurse_decks(draw(d, c), total + c);
    }
    return id;
}

//precalculate all possible deck state for first three cards before main calculation begins
static void all_subdecks(const Deck& d) 
{
    for (int i = 1; i < 11; i++) 
    {
        if (d[i] == 0) continue;
        Deck d1 = draw(d, i);
        for (int j = 1; j < 11; j++) 
        {
            if (d1[j] == 0) continue;
            Deck d2 = draw(d1, j);
            for (int k = 1; k < 11; k++) 
            {
                if (d2[k] == 0) continue;
                Deck rem = draw(d2, k);
                recurse_decks(rem, i + j);
                if (i == j) recurse_decks(rem, i); 
            }
        }
    }
}

static inline void print_line_value(const char* label, double v) 
{
    if (!std::isfinite(v)) printf("%s: null\n", label);
    else                   printf("%s: %.12f\n", label, v);
}

//iterate through every possible starting the player could be dealt and every possible dealer card, then calculate total EV of entire shoe
static void compute_and_print_ev(const Deck& deck_rem) {
    double total_ev = 0.0;
    double insurance_ev_total = 0.0;
    double insurance_offer_p  = 0.0;

    //pre-map all reachable rem-deck states and needed split depth.
    all_subdecks(deck_rem);

    #pragma omp parallel for reduction(+:total_ev,insurance_ev_total,insurance_offer_p)
    for (int i = 1; i < 11; i++) 
    {
        if (deck_rem[i] <= 0) continue;
        double p1 = (double)deck_rem[i] / deck_rem[0];
        Deck d1 = draw(deck_rem, i);

        for (int j = i; j < 11; j++) 
        {
            if (d1[j] <= 0) continue;
            double p2 = (double)d1[j] / d1[0];
            if (i != j) p2 *= 2.0;
            Deck d2 = draw(d1, j);

            for (int k = 1; k < 11; k++) {
                if (d2[k] <= 0) continue;
                double p3 = (double)d2[k] / d2[0];
                Deck rem = draw(d2, k);

                int d_idx = 0;
                auto it = indexer.find(rem);
                if (it != indexer.end()) d_idx = it->second;
                else 
                {
                    continue;
                }

                ActionEVs a = compute_action_evs(d_idx, i, j, k);
                double best = a.player_blackjack ? a.stand : best_ev_value(a);

                const double deal_p = p1 * p2 * p3;
                total_ev += deal_p * best;

                if (k == 1) 
                {
                    insurance_offer_p += deal_p;
                    const double p_hole_ten = (rem[0] > 0) ? ((double)rem[10] / rem[0]) : 0.0;
                    const double insurance_ev_per_ins_bet = 3.0 * p_hole_ten - 1.0;
                    insurance_ev_total += deal_p * insurance_ev_per_ins_bet;
                }
            }
        }
    }

    double player_ev_pct = total_ev * 100.0;
    double ins_cond_pct = (insurance_offer_p > 0.0) ? (insurance_ev_total / insurance_offer_p) * 100.0 : 0.0;

    print_line_value("player_ev_pct", player_ev_pct);
    print_line_value("insurance_ev_pct", ins_cond_pct);
    fflush(stdout);
}

//print EV for every action for a specific dealer up card and player hand
static void compute_and_print_hand(const Deck& deck_rem, int up, int p1, int p2) {
    // deck_rem is the remaining shoe AFTER removing visible cards (up, p1, p2).
    require_card_val(up, "dealer upcard");
    require_card_val(p1, "player card 1");
    require_card_val(p2, "player card 2");

    const int base_total = (p1 == p2) ? p1 : (p1 + p2);
    int d_idx = recurse_decks(deck_rem, base_total);

    ActionEVs a = compute_action_evs(d_idx, p1, p2, up);

    double best_ev = NAN;
    const char* best_nm = a.player_blackjack ? "Stand" : best_action_name(a, best_ev);
    if (a.player_blackjack) best_ev = a.stand;

    auto ev_or_null = [](double v) -> std::string 
    {
        if (!std::isfinite(v)) return "null";
        char buf[64];
        snprintf(buf, sizeof(buf), "%.12f", v);
        return std::string(buf);
    };

    printf("up=%d p1=%d p2=%d\n", up, p1, p2);
    printf("stand: %s\n", ev_or_null(a.stand).c_str());
    printf("hit: %s\n", ev_or_null(a.hit).c_str());
    printf("double: %s\n", ev_or_null(a.dbl).c_str());
    printf("split: %s\n", ev_or_null(a.split).c_str());
    printf("best: %s, %s\n", best_nm, ev_or_null(best_ev).c_str());
    fflush(stdout);
}

static int parse_int(const char* s, const char* label) 
{
    char* end = nullptr;
    long v = std::strtol(s, &end, 10);
    if (!s || *s == '\0' || !end || *end != '\0') fatalf("Bad %s: %s", label, s ? s : "(null)");
    if (v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max())
        fatalf("Out of range %s: %s", label, s ? s : "(null)");
    return (int)v;
}

static void remove_visible_cards(Deck& deck, int up, int p1, int p2) 
{
    const int cards[3] = {up, p1, p2};
    for (int i = 0; i < 3; i++) {
        int c = cards[i];
        require_card_val(c, "visible card");
        if (deck[c] <= 0) fatalf("Not enough cards for value %d", c);
        deck[c] -= 1;
        deck[0] -= 1;
    }
}

int main(int argc, char** argv) 
{
    initialize();
    dealer_sequences();

    if (argc < 2) 
    {
        fatalf("Usage: blackjack EV <c1..c10> | blackjack HAND <c1..c10> <up> <p1> <p2>");
    }

    std::string cmd = argv[1];

    if (cmd == "EV") 
    {
        if (argc != 12) fatalf("EV requires 10 counts: EV <c1..c10>");
        Deck deck = {0};
        for (int c = 1; c <= 10; c++) 
        {
            int v = parse_int(argv[1 + c], "count");
            deck[c] = v;
            deck[0] += v;
        }
        compute_and_print_ev(deck);
        return 0;
    }

    if (cmd == "HAND") 
    {
        if (argc != 15) fatalf("HAND requires 10 counts + 3 cards: HAND <c1..c10> <up> <p1> <p2>");
        Deck deck = {0};
        for (int c = 1; c <= 10; c++) 
        {
            int v = parse_int(argv[1 + c], "count");
            deck[c] = v;
            deck[0] += v;
        }
        int up = parse_int(argv[12], "dealer upcard");
        int p1 = parse_int(argv[13], "player card 1");
        int p2 = parse_int(argv[14], "player card 2");

        remove_visible_cards(deck, up, p1, p2);
        compute_and_print_hand(deck, up, p1, p2);
        return 0;
    }

    fatalf("Unknown command: %s", cmd.c_str());
}
