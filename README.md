CD Calculator For Pragmatic Play live blackjack table rules:

```
8 deck, S17, DAS, 1 resplit only, Hybrid ENHC (Dealer checks for BJ on ace, but not 10)
```

To compile:

```
g++ -O3 -fopenmp blackjackCDC.cpp -o blackjackCDC
```

To calculate the EV of a specific deck composition
```
blackjackCDC EV <deck composition>
```

