
# markov chain of 7 states : posible genopytpe structure of borher and sister
# 4abcd - 4 alleles of a locus, both individuals have two different alleles
# 3aabc, 3abac, - 3 alleles of a locus (either homozyg and heterozyg, or 2 heterozyg)
# 2abab, 2aaab, 2aabb, - 2 alleles (2x heter, 1 homo 1 hetero, 2 homo)
# 1aaaa - homogenised locus
mc <- as.matrix(read.table('data/inbreeding_mc.tsv', sep = '\t'))

# sanity check
# all(rowSums(mc) == 1)
# [1] TRUE

# absorb_states is the state that lock MC in that state : i.e. transition probability to itself is 1
absorb_states <- which(diag(mc) == 1)

# number of states (n) and number of non absorb states (m)
n <- nrow(mc)
m <- n - length(absorb_states)

# subsampling matrix to its only non-absobing states
non_abs_mc <- mc[-absorb_states, -absorb_states]

# calculating fundamental matrix as inversion of (E - non_abs_mc), where E is diagonal unit matrix of same dimension as non_abs_mc
fundamental_matrix <- solve(diag(rep(1, m)) - non_abs_mc)

# calculate mean number of steps to get into absorbing state given initial state
fundamental_matrix %*% rep(1, m)
# 7.666667 6.666667 7.166667 5.666667 4.833333 6.666667

#######
# I am interested in probabilities of getting into absorbing states for every timepoint and initial state
# not just the mean number of steps
# to do that, I calculate it using brure force - a while loop
# following function calculates distribution of probability to be in absorbing state given initial state

calc_distr <- function(mc, t0, treshold = 0.999){
    absorb_state <- which(diag(mc) == 1)
    state_probabilities <- t0 %*% mc
    psr_of_being_absorbed <- state_probabilities[absorb_state]
    psr <- 0
    while(psr < treshold){
        state_probabilities <- state_probabilities %*% mc
        psr <- state_probabilities[absorb_state]
        psr_of_being_absorbed <- c(psr_of_being_absorbed, psr)
    }
    return(psr_of_being_absorbed)
}

# complete distribution for m inital states
# distributions[[1]] ... AB  CD
# distributions[[2]] ... AB  AC
# distributions[[3]] ... AA  BC
# distributions[[4]] ... AB  AB
# distributions[[5]] ... AA  AB
# distributions[[6]] ... AA  BB
# 7th state : AA AA is the absorbing state

distributions <- list()
for (i in 1:m) {
    t0 <- rep(0, n)
    t0[i] <- 1
    distributions[[i]] <- calc_distr(mc, t0)
}

png('figures/heterozigosity_loss_during_inbreeding.png')

plot(NULL,
     xlim = c(0, max(unlist(lapply(distributions, length)))),
     ylim = c(0, 1),
     xlab = "inbreeding generation",
     ylab = 'probability of homozygous state',
     main = 'loss of heterozygosity during inbreeding')

palette <- rainbow(m, s = 0.5)
for (i in 1:m){
    lines(distributions[[i]], lwd = 1.8,
          col = palette[i])
    points(distributions[[i]], pch = 20, cex = 0.8,
           col = palette[i])
}
grid()

legend('bottomright',
       legend = c('AB  CD', 'AB  AC', 'AA  BC', 'AB  AB', 'AA  AB', 'AA  BB'),
       pch = 20, col = palette, title = 'initial allelic state', bty = 'n')


dev.off()

##### Recalculate worst case scenario with various numbers of independent loci
# a, all loci start at AB CD

png('figures/heterozigosity_loss_num_of_ABCD_loci.png')

plot(NULL,
     xlim = c(0, max(unlist(lapply(distributions, length)))),
     ylim = c(0, 1),
     xlab = "inbreeding generation",
     ylab = 'probability of fully homozygous loci',
     main = 'heterozygosity loss relative to number of indipendent AB CD loci')

number_of_loci <- c(1,5,10,20,100)
palette <- rev(heat.colors(length(number_of_loci)))
for (i in 1:length(number_of_loci)){
    proportions <- distributions[[1]]^number_of_loci[i]
    lines(proportions, lwd = 1.8, col = palette[i])
    points(proportions, pch = 20, cex = 0.8, col = palette[i])
}
grid()

legend('bottomright',
       legend = number_of_loci,
       pch = 20, col = palette,
       title = 'number of indipendent loci', bty = 'n')

dev.off()

# b, all starts at AA AB

png('figures/heterozigosity_loss_num_of_AAAB_loci.png')

plot(NULL,
     xlim = c(0, max(unlist(lapply(distributions, length)))),
     ylim = c(0, 1),
     xlab = "inbreeding generation",
     ylab = 'probability of fully homozygous loci',
     main = 'heterozygosity loss relative to number of indipendent AA AB loci')

number_of_loci <- c(1,5,10,20,100)
palette <- rev(heat.colors(length(number_of_loci)))
for (i in 1:length(number_of_loci)){
    proportions <- distributions[[5]]^number_of_loci[i]
    lines(proportions, lwd = 1.8, col = palette[i])
    points(proportions, pch = 20, cex = 0.8, col = palette[i])
}
grid()

legend('bottomright',
       legend = number_of_loci,
       pch = 20, col = palette,
       title = 'number of indipendent loci', bty = 'n')

dev.off()