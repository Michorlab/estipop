
# Specify two fixed transitions, birth and death
transitionList = TransitionList(2, FixedTransition(population = 0, rate = Rate(params[1]*exp(-t)), fixed = c(2)),
                                FixedTransition(population = 0, rate = Rate(params[2]*exp(-t)), fixed = c(0)))

estimateBP(seq(0,5),1,transitionList, matrix(c(1,2,3,4,5),ncol=1), c(1,1))
