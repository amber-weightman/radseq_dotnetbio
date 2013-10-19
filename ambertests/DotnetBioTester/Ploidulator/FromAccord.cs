using Accord.Statistics.Distributions.Univariate;
using Accord.Statistics.Models.Markov;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Globalization;

namespace Ploidulator
{
    // num choc chips in a biscuit is measured at 100
    // what is the probability that the first half of the biscuit will contain 80% of the choc chips?




    // num events A required to get num items B

    // for each cluster the mean is
    // mean sequence frequency (where fx is the num distinct individuals which have that seq at least once)
    // where frequency is

    // for each indivdual in the cluster (make the individual a sub cluster)
    // for all the sequences attributed to that individual, what is the mean seq frequency
    //[553, 334, 120, 52, 5, 2, 1, 1, 1]

    // WHAT IS THE PROBABILITY THAT AN INDIVIDUAL HAS EXACTLY TWO DISTINCT SEQUENCES?
    // WHERE ON AVERAGE FOR THIS CLUSTER, EACH INDIVIDUAL HAS X DISTINCT SEQUENCES
    // -- but I still don't get it. We can count all these values. Where is the unknown?

    // rate the individual based on how accurately its sequences represent top 2 for itself
    // rate the individual based on how accurately its top 2 sequences mirror top 2 for the group
    // if the individual rates low against the group, we can exclude the individual?
    // by removing the individual this will improve the group rating



    ///////////////////////////////////////////////////////////////////////////////////////////////
    // A. setup
    // ploidy = N

    // 1. compute for cluster...
    // A cluster has M sequences
    // A cluster has N top sequences
    // The closer N is to M the better the cluster is
    // But also, if N!= M, N+1 should be significantly small

    // 2. Compute for first individual, then if necessary recompute for cluster...
    // Each individual has P sequences
    // Each individual has <= N top sequences
    // If the <=N top sequences != the cluster's N sequences, the individual does not belong. Remove it
    // Otherwise, the closer N is to M the better the individual's ploidy is

    // Or, rather than removing from cluster, call it a degree of copy number variation

    // 3. Repeat for all remaining individuals

    // This would not exclude individuals based on minor sequences, but it would exclude individuals which
    // do not match the genotype of the group

    // Q: Can we exclude individuals which do not represent the ploidy of the cluster they have been
    // placed in (exclude all reads from that individual from that cluster)? This will still retain the
    // outlying sequences.


    // Or, compute the accuracy for each individual (% reads within ploidy)
    // The cluster accuracy is the combined accuracy of its individuals

    ////////////////////////////////////////////////////////////////////////////////////////////////




    // mean: add up and divide by total (λ)
    // median: middle value (or the mean of the two middle values) [5]
    // mode: most frequently repeated value [1]
    // range: diff between max and min

    class FromAccord
    {
        //The Poisson distribution is a discrete probability distribution that expresses the 
        //probability of a number of events occurring in a fixed period of time if these events 
        //occur with a known average rate and independently of the time since the last event.

        // poisson distribution is constructed using lambda
        // k is used to calculate distribution or mass - it can be passed in for various values
        // p can also be passed in for various values
        // x can also be passed in for various values

        // one of these to be created for each cluster
        // what is lambda?

        // event occurrs λ times over a specified time interval
        // what is the probability of exactly x num occurrences during that period?

        // e.g. normally we get 10 customers over time[t], what is the probability
        // of 7 customers over time[t]?

        // λ is the number of occurrences in the time frame we are talking about
        // mu is another name for λ (woopsy u
        // variance == λ (o with hat
        // mean == λ
        // mu == mean == expected value


        // x must be whole num
        // x is a random var representing the num occurrences
        // p(x=1) .. prob that x is 1

        // only input is λ and x

        // num choc chips in a cookie


        // binomial   n = size of sample     p = probability of 1 occurrence 
        // whereas λ = n*p
        // large n and small p makes poisson reasonable   n>50 && np < 5 very roughly

        // if n and p are unknown, poisson it. we only need to know the mean to use poisson
        // although we should also know that n is v large and p is v small
        // (poisson approximates binomial)

        // the reason poisson is slightly inaccurate is that the mean is not always equal 
        // to the variance

        // negative binomial? generalised poisson?

        // i think hmm was multinomial?

        public void PoissonDist(double lambdaVal)
        {
            // Create a Poisson distribution with λ = 4.2 
            PoissonDistribution dist = new PoissonDistribution(lambda: lambdaVal);

            // Common measures 
            double mean = dist.Mean;     // 4.2 
            double median = dist.Median; // 4.0 
            double var = dist.Variance;  // 4.2 

            // Cumulative distribution functions 
            double cdf = dist.DistributionFunction(k: 2);               // 0.39488100648845126 
            double ccdf = dist.ComplementaryDistributionFunction(k: 2); // 0.60511899351154874 

            // Probability mass functions 
            double pmf1 = dist.ProbabilityMassFunction(k: 4); // 0.19442365170822165 
            double pmf2 = dist.ProbabilityMassFunction(k: 5); // 0.1633158674349062 
            double pmf3 = dist.ProbabilityMassFunction(k: 6); // 0.11432110720443435 
            double lpmf = dist.LogProbabilityMassFunction(k: 2); // -2.0229781299813 

            // Quantile function 
            int icdf1 = dist.InverseDistributionFunction(p: 0.17); // 2 
            int icdf2 = dist.InverseDistributionFunction(p: 0.46); // 4 
            int icdf3 = dist.InverseDistributionFunction(p: 0.87); // 7 

            // Hazard (failure rate) functions 
            double hf = dist.HazardFunction(x: 4); // 0.19780423301883465 
            double chf = dist.CumulativeHazardFunction(x: 4); // 0.017238269667812049 

            // String representation 
            string str = dist.ToString(System.Globalization.CultureInfo.InvariantCulture); // "Poisson(x; λ = 4.2)"

        }


        //The negative binomial distribution is a discrete probability distribution of the number 
        //of successes in a sequence of Bernoulli trials before a specified (non-random) number of 
        //failures (denoted r) occur. For example, if one throws a die repeatedly until the third 
        //time “1” appears, then the probability distribution of the number of non-“1”s that had 
        //appeared will be negative binomial.

        // where a success is a ploidy of 2 and a failure is anything else?
        // neg binomial is a generalisation of geometric distribution
        // (num coin tosses required to get the r-th head)
        // prob that the Ath event occurrs on the Bth trial

        // binomial - get the num successes x ina fixed number of trials n
        // neg binomial - get the number of trials (number of failures) x needed to get a fixed number of successes r
        // where for both above, x is the random var
        // P success = p (l.c.)
        // P failure = 1-p
        // X (u.c.) represents trial num of the r-th success

        private void NegBinomial()
        {
            // Create a new Negative Binomial distribution with r = 7 and p = 0.42 
            var dist = new NegativeBinomialDistribution(failures: 7, probability: 0.42);

            // Common measures 
            double mean = dist.Mean;     // 5.068965517241379 
            double median = dist.Median; // 5.0 
            double var = dist.Variance;  // 8.7395957193816862 

            // Cumulative distribution functions 
            double cdf = dist.DistributionFunction(k: 2);               // 0.19605133020527743 
            double ccdf = dist.ComplementaryDistributionFunction(k: 2); // 0.80394866979472257 

            // Probability mass functions 
            double pmf1 = dist.ProbabilityMassFunction(k: 4); // 0.054786846293416853 
            double pmf2 = dist.ProbabilityMassFunction(k: 5); // 0.069908015870399909 
            double pmf3 = dist.ProbabilityMassFunction(k: 6); // 0.0810932984096639 
            double lpmf = dist.LogProbabilityMassFunction(k: 2); // -2.3927801721315989 

            // Quantile function 
            int icdf1 = dist.InverseDistributionFunction(p: 0.17); // 2 
            int icdf2 = dist.InverseDistributionFunction(p: 0.46); // 4 
            int icdf3 = dist.InverseDistributionFunction(p: 0.87); // 8 

            // Hazard (failure rate) functions 
            double hf = dist.HazardFunction(x: 4); // 0.10490438293398294 
            double chf = dist.CumulativeHazardFunction(x: 4); // 0.64959916255036043 

            // String representation 
            string str = dist.ToString(CultureInfo.InvariantCulture); // "NegativeBinomial(x; r = 7, p = 0.42)"
        }


        //The arbitrary-density Hidden Markov Model uses any probability density function 
        //(such as GaussianMixture Model) for computing the state probability. In other words, 
        //in a continuous HMM the matrix of emission probabilities B is replaced by an array of 
        //either discrete or continuous probability density functions.
        private void Hmm()
        {
            // Create the transation matrix A 
            double[,] transitions = { { 0.7, 0.3 }, { 0.4, 0.6 } };

            // Create the vector of emission densities B
            GeneralDiscreteDistribution[] emissions = 
            {  
                new GeneralDiscreteDistribution(0.1, 0.4, 0.5),
                new GeneralDiscreteDistribution(0.6, 0.3, 0.1)
            };

            // Create the initial probabilities pi 
            double[] initial = { 0.6, 0.4 };

            // Create a new hidden Markov model with discrete probabilities 
            var hmm = new HiddenMarkovModel<GeneralDiscreteDistribution>(transitions, emissions, initial);

            // After that, one could, for example, query the probability 
            // of a sequence ocurring. We will consider the sequence 
            double[] sequence = new double[] { 0, 1, 2 };

            // And now we will evaluate its likelihood 
            double logLikelihood = hmm.Evaluate(sequence);

            // At this point, the log-likelihood of the sequence 
            // ocurring within the model is -3.3928721329161653. 

            // We can also get the Viterbi path of the sequence 
            int[] path = hmm.Decode(sequence, out logLikelihood);

            // At this point, the state path will be 1-0-0 and the 
            // log-likelihood will be -4.3095199438871337        
        }


        public void CreateGMM()
        {
            // normal == gaussian distribution
            NormalDistribution[] distributions = new NormalDistribution[] 
            {
                new NormalDistribution(2, 1), new NormalDistribution(5, 1)
            };

            //MixtureModel(distributions);
        }

        public void CreatePMM()
        {
            // normal == gaussian distribution
            PoissonDistribution[] distributions = new PoissonDistribution[] 
            {

                new PoissonDistribution(lambda: 4.2), new PoissonDistribution(lambda: 4.5)
            };

            double[] pmfs = MixtureModel(distributions);
        }

        // Create a new mixture containing Poisson distributions
        public Mixture<PoissonDistribution> GetMixture(PoissonDistribution[] distributions)
        {
            return new Mixture<PoissonDistribution>(distributions);
        }
        // Create a new mixture containing Normal (Gaussian) distributions
        public Mixture<NormalDistribution> GetMixture(NormalDistribution[] distributions)
        {
            return new Mixture<NormalDistribution>(distributions);
        }



        public double[] MixtureModel(PoissonDistribution[] distributions)
        {
            Mixture<PoissonDistribution> mix = new Mixture<PoissonDistribution>(distributions);

            // Common measures 
            double mean = mix.Mean;     // 3.5 
            double median = mix.Median;   // 3.4999998506015895 
            double var = mix.Variance; // 3.25 

            // Cumulative distribution functions 
            double cdf = mix.DistributionFunction(x: 4.2);               // 0.59897597553494908 
            double ccdf = mix.ComplementaryDistributionFunction(x: 4.2); // 0.40102402446505092 

            // Probability mass functions 
            double pmf1 = mix.ProbabilityDensityFunction(x: 1.2); // 0.14499174984363708 
            double pmf2 = mix.ProbabilityDensityFunction(x: 2.3); // 0.19590437513747333 
            double pmf3 = mix.ProbabilityDensityFunction(x: 3.7); // 0.13270883471234715 
            double lpmf = mix.LogProbabilityDensityFunction(x: 4.2); // -1.8165661905848629 

            // Quantile function 
            double icdf1 = mix.InverseDistributionFunction(p: 0.17); // 1.5866611690305095 
            double icdf2 = mix.InverseDistributionFunction(p: 0.46); // 3.1968506765456883 
            double icdf3 = mix.InverseDistributionFunction(p: 0.87); // 5.6437596300843076 

            // Hazard (failure rate) functions 
            double hf = mix.HazardFunction(x: 4.2);            // 0.40541978256972522 
            double chf = mix.CumulativeHazardFunction(x: 4.2); // 0.91373394208601633 

            // String representation: 
            // Mixture(x; 0.5 * N(x; μ = 5, σ² = 1) + 0.5 * N(x; μ = 5, σ² = 1)) 
            string str = mix.ToString(CultureInfo.InvariantCulture);

            return new double[] { pmf1, pmf2, pmf3 };
        }

    }
}
