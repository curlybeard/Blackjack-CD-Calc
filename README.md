**Composition Dependent Calculator For Pragmatic Play live blackjack table rules:**

```
8 deck, S17, DAS, 1 resplit only, Hybrid ENHC (Dealer checks for BJ on ace, but not 10)
```

**To compile:**

```
g++ -O3 -fopenmp blackjackCDC.cpp -o blackjackCDC
```

**USAGE:**

**1) To calculate the EV of a specific deck composition:**

```
blackjackCDC EV <deck composition>
```

eg) Calculate EV of a full 8 deck shoe
```
blackjackCDC EV 32 32 32 32 32 32 32 32 32 128
```

**2) To calculate the EV of all player actions given a player hand and dealer up card:**

```
blackjackCDC HAND <deck composition> <dealer up card> <player card 1> <player card 2>
```

eg) Calculate the optimal decision of a player hand of a 10 and 6, versus a dealer 10, starting from a full 8 deck shoe
```
blackjackCDC HAND 32 32 32 32 32 32 32 32 32 128 10 10 6
```

**TIME COMPLEXITY = O(n)** 

O(V * D), where D is the number of unique deck states, and V is the number of card face values (1-10)
