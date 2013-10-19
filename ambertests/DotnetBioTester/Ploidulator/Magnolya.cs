/*using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Bio.Algorithms.Metric
{

    private enum Gamma {None, Haploid, Diploid, Triploid, Uniform};

    // "unique sequence in cluster"
    private class MagnolyaContig
        {
            private string contigId;    // being the query string
            private int contigLength;   // query string length
            private int countSample1;   // num samples for indiv
            private int countSample2;   // num samples for indiv
            
            //private string clusterId;
            private int start;          // start index
            private int end;            // end index
            //private int clusterLen;

            public string ContigId { get; set; }
            public int ContigLength { get; set; }
            public int CountSample1 { get; set; }
            public int CountSample2 { get; set; }
            public int Start { get; set; }
            public int End { get; set; }
        }


    private class PoissonMixture
    {
        #region vars
                private string clusterId;   // cluster/chromosome id
                public string ClusterId { get; set; }
                private List<MagnolyaContig> contigs;
                public List<MagnolyaContig> Contigs { get; set; }

                // "SELF"
                private List<int> mu;public List<int> Mu {get; set;}
                private List<List<int>> mus;public List<List<int>> Mus {get; set;}
                private List<int> mix;public List<int> Mix {get; set;}
                private List<double> bics;public List<double> BICs {get; set;}
                private List<List<int>> mixs;public List<List<int>> Mixs {get; set;}
                private int l;public int L {get; set;}
                private List<int> likelihoods;public List<int> Likelihoods {get; set;}

                private int alpha;public int Alpha {get; set;}
                private List<int> alphas;public List<int> Alphas {get; set;}

                private List<double> bics;public List<double> Bics {get; set;}
                private List<float> x;public List<float> X {get; set;}

        
                // constructor must set contigs
        #endregion
        private Gamma ploidy;
        private List<int> counts;
        private List<int> intervals;
        List<int> k;
        double epsilon;
        int max_iters;



        public PoissonMixture(List<int> counts, List<int> intervals, List<int> k = null, 
            double epsilon = 1e-15, int max_iters = 100, Gamma gamma = Gamma.None)
        {
            this.ploidy = gamma;
            this.counts = counts;
            this.intervals = intervals;
            this.k = k;
            this.epsilon = epsilon;
            this.max_iters = max_iters;
            Console.WriteLine("Gamma = "+gamma);
        
            EmPoisson(counts, intervals,k, epsilon,max_iters,gamma); // sets L, mu, mix

            Console.WriteLine("L: " + this.L);
            Console.WriteLine("lambda: " + this.mu);
            Console.WriteLine("p(k): " + this.mix);
        }


    
        private void SetModelNumber(int KsIndex)
        {
            L = Likelihoods[KsIndex];
            Mu = Mus[KsIndex];
            Mix = Mixs[KsIndex];
            Alpha = Alphas[KsIndex];

            Console.WriteLine("Changed model parameters to:\nL: "+L.ToString() 
                + "\nlambdas: "+Mu + "\np(k): "+Mix);
        }
        
        private float PriorDepth(int x, int t)
        {
            float px = 0;
            for(int k = 0; k < Mix.Count; k++)
            {
                px+=Mix[k] * Pmf(x, Mu[k], t);
            }
            return px;
        }


        //Posterior probability:
        //p(k|x) = p(k). (px|k).p(x))
        private List<float> Posterior(int x, int t)
        {
            List<float> posterior = new List<float>();
            float px = PriorDepth(x,t);

            for(int k = 0; k < Mix.Count; k++)
            {
                posterior.Add( Mix[k] * Pmf(x, Mu[k], t) / px);
            }
            return posterior;
        }
        
        //evaluate gaussian at x
        private int Pmf(int x, int mu, int t)
        {
            //scripy
            return stats.poisson.pmf(x, mu*t);
        }

        // instead of returnign sets values of l, mu and mix directly
        //THERE SHOULD BE 2 OF THESE METHODS
        public void EmPoisson(List<int> counts, List<int> intervals, List<int> k, double epsilon, 
            int maxIters, Gamma gamma)
        {
            // remove contigs with 0 length
            //counts    = array([counts[i]    for i in range(len(intervals)) if intervals[i]!=0])
            //intervals = array([intervals[i] for i in range(len(intervals)) if intervals[i]!=0])
            foreach(int num in intervals)
            {
                if(num != 0)
                {
                    counts.Add(num);
                }
            }
            intervals = counts;

            // filter repeats with too high copy number
            X = H.Divide(counts, intervals);

            Console.WriteLine("k: "+k);
    
            // sets [self.likelihoods, self.mus, self.mixs, self.BICs, self.ks, self.alphas] 
            EmPoissonBic(counts, intervals,k,epsilon, max_iters);
                              
            int best = Bics.IndexOf(Bics.Min());
            this.L = Likelihoods[best];
            this.Mu = Mus[best]; 
            this.Mix = Mixs[best];
        }


        //// sets [self.likelihoods, self.mus, self.mixs, self.BICs, self.ks, self.alphas] 
                // k = range of cluster numbers in a list
        public void EmPoissonBic(List<int> count, List<int> interval, 
                List<int> k = null, double epsilon = 2e-13, 
                int max_iters = 100, Gamma gamma = Gamma.Diploid, string figname = "poisson.pdf", int addUniform = 1,
                List<double> initPercLists = null, bool bw = true, int plotContigSize = 1000, int maxVal = int.MaxValue, 
                List<int> ymax = null)
        {
            if(count.Count == 0) { Console.WriteLine("Warning: Zero reads detected"); }

            Console.WriteLine("****Starting model selection using BIC****");
            Console.WriteLine("Number of contigs used to train EM: "+count.Count);
            Console.WriteLine("Median contig length: ", H.Median(interval));
    
            k = k == null ? new List<int>(new int[]{3,4,6,8,10,12,15,20}) : k;
            
            List<int> likelihoods = new List<int>();    List<int> mus = new List<int>();
            List<int> mixs = new List<int>();           List<int> alphas = new List<int>();

            foreach(int sequence in k) // for each cluster
            {
                Dictionary<string, int> results = EmPoissonInit(count, interval, 
                        sequence, 0, epsilon, max_iters, gamma, 
                        figname, addUniform, initPercLists, plotContigSize, maxVal);
                likelihoods.Add(results["likelihood"]);
                mus.Add(results["mu"]);
                mixs.Add(results["mix"]);
                alphas.Add(results["alpha"]);
            }

            // Calculate Bayesian Information Criterion
            int n = count.Count;
            Bics = new List<double>();

            for(int i = 0; i < likelihoods.Count; i++)
            {
                double val = (-2 * likelihoods[i] + (k[i]+1+2*addUniform) * Math.Log(n)); // pi's + lambda + geom (alpha and pi);
                Bics.Add(val);
                Console.WriteLine(String.Format("k = %d; BIC = %2.2f; L = %2.2f; lambda = %2.2e;",
                          k[i], Bics[i], likelihoods[i], mus[i])); // was mus[i][0]    
            }
            int best = Bics.IndexOf(Bics.Min());
            this.Likelihoods = likelihoods;
            this.Mus = mus;
            this.Mixs = mixs; // todo fixme I need to resolve this
            this.k = k;
            this.alphas = alphas;

            this.L = Likelihoods[best];
            this.Mu = Mus[best]; 
            this.Mix = Mixs[best];
        }
    

        public Dictionary<string, int> EmPoissonInit(List<int> count, List<int> interval, 
                int k = 10, int plot_mixture = 0, double epsilon = 2e-13, int max_iters = 100, 
                Gamma gamma = Gamma.Haploid, string figname = "poisson.pdf", int addUniform = 1, 
                List<double> initPercList = null, int plotContigSize = 1000, int maxval = int.MaxValue)
        {
            //INITPERCLIST    List with percentiles to initialize the EM
            if(initPercList == null && (gamma != Gamma.Haploid && gamma != Gamma.Diploid))
            {
                initPercList = new List<double>(new double[]{1, 5, 12.5, 25, 50, 70});   
                Console.WriteLine("Performing initializations at 1, 5, 12.5, 25 and 50th percentile");
            } 
            else if(initPercList == null && (gamma == Gamma.Haploid || gamma == Gamma.Diploid))
            {
                initPercList = new List<double>(new double[]{50});   
                Console.WriteLine("Performing initialization at 50th percentile");
            }
            
            List<int> likelihoods = new List<int>();    List<int> mus = new List<int>();
            List<int> mixs = new List<int>();           List<int> alphas = new List<int>();

            foreach(double perc in initPercList) 
            {
                // lambda for 1, 5, 12.5, 25, 50, 70 is calculated below (currently returns null though)
                List<double> lambda_init = ScoreAtPercentile(H.Divide(count, interval), perc); // todo fixme


                List<int> results = EmPoissonB(count,interval,k,0,epsilon,max_iters,
                        gamma, figname, addUniform, lambda_init, plotContigSize, maxval);
                int likelihood = results[0];
                int mu = results[1];
                int mix = results[2];
                int alpha = results[3];

                likelihoods.Add(likelihood);
                mus.Add(mu);
                mixs.Add(mix);
                alphas.Add(alpha); // these are all local copies, not class copies
                
                Console.WriteLine(String.Format("Lambda init: %2.2e for %2.2f-th percentile -> L = %2.2f  ",
                          lambda_init, perc,L)); 
            }
            int best = likelihoods.IndexOf(likelihoods.Max());
            
            Dictionary<string, int> returnDict = new Dictionary<string, int>();   
            returnDict.Add("likelihood", likelihoods[best]);
            returnDict.Add("mu", mus[best]);
            returnDict.Add("mix", mixs[best]);
            returnDict.Add("alpha", alphas[best]);
            return returnDict;
        }
    

        
         //[L, MU, S2, MIX] = EM (COUNT, INTERVAL, K, STEPK, FIX_S2, REG, PLOT, EPS, MAXITER)
    
    //Uses the EM algorithm to fit a mixture of K 1D Poissons, where the means
    //and the variances of the models are spaced as (1 + k/STEPK), k = 0..K-1.
    //Optional parameters are: 
    //FIX_S2	set the variance to the mean (1, dflt.) or estimate (0)
    //REG     regularisation constant (dflt. 0)
    //PLOT    plot progress (dflt. 0)
    //EPS     epsilon, the change in likelihood to stop training (dflt. 1e-5)
    //MAXITER maximum number of iterations (dlft. 100)
    //Returns the likelihood L, the means MU(), the variances S2() and priors MIX.
         
        public List<int> EmPoissonB(List<int> count, List<int> interval, int k = 10, int plot_mixture = 0, 
            double epsilon = 2e-13, int max_iters = 100, Gamma gamma = Gamma.Haploid, string figname = "poisson.pdf", 
            int addUniform = 1,
            List<double> lambda_init = null, int plotContigSize = 1000, int maxval = int.MaxValue)
        {

            if(addUniform != 0 && addUniform != 1) { throw new Exception("addUniform should be 0 or 1"); }

            
            k    = k+addUniform;  // One additional k for the Uniform distribution 
            int[] x    = count.ToArray();

            int[] t    = interval.ToArray();
            float N = (float)x.Length;
            float A = (float)k;
            // these should prob just be arrays
            List<double> mu = new List<double>(new double[k]); // all to have value 0.0 assigned
            List<double> mix = new List<double>(new double[k]); // all to have value 0.0 assigned
            List<double> sumR = new List<double>(new double[k]); // all to have value 0.0 assigned
            double alpha = 0.0001;
            List<double> mix_new = new List<double>(new double[k]); // all to have value 0.0 assigned
            List<double> mu_new = new List<double>(new double[k]); // all to have value 0.0 assigned
                        
            double min_data = (H.Divide(x,t).Min()) > epsilon ? H.Divide(x,t).Min() : epsilon;
            float max_data = H.Divide(x,t).Max(); //# Normalised for interval
            int tries = 0; 
            int retry = 1;
            float[] gammaArray;


            while(retry == 1 && tries <= 5)
            {
                // todo fixme i don't think these are all right
                if(gamma == Gamma.Haploid || gamma == Gamma.Diploid) // haploid
                {
                    gammaArray = Dirichlet(gamma,k,N);
                    if(lambda_init == null)
                    {
                        lambda_init = new List<double>(); // should just be double
                        lambda_init.Add(H.Median(H.Divide(x,t)) / (int)gamma);
                        Console.WriteLine(String.Format("Haploid/Diploid initialization, lamda: " , lambda_init.ToString()));
                    }
                    for(int i = 0; i < A; i++)
                    {
                        mu[i]  = H.Multiply(lambda_init, (i+1))[0];
                        mix[i] = 0.4/(k-1);
                    }            
                    mix[0] = 0.6;
                }
                else if(gamma == Gamma.Triploid) // allo-triploid
                {
                    gammaArray = Dirichlet(gamma,k,N);
                    if(lambda_init == null)
                    {
                        lambda_init = new List<double>();
                        lambda_init.Add(H.Median(H.Divide(x,t)) / (int)gamma);
                    }
                    for(int i = 0; i < A; i++)
                    {
                        mu[i]  = H.Multiply(lambda_init, (i+1))[0];
                        mix[i] = 0.2/(k-2);
                    }
                    mix[1] = 0.4;
                    mix[2] = 0.4;   
                    Console.WriteLine("Triploid initialization, lamda: " + lambda_init.ToString());
                }
                else if(gamma == Gamma.Uniform) // uniform
                {
                    gammaArray = Dirichlet(gamma,k,N);
                    if(lambda_init == null)
                    {
                        lambda_init = new List<double>();
                        lambda_init.Add(H.Median(H.Divide(x,t)));
                        Console.WriteLine("Default initialization, lamda: " + lambda_init.ToString());
                    }
                    for(int i = 0; i < A; i++)
                    {
                        mu[i]  = H.Multiply(lambda_init, (i+1))[0];
                        mix[i] = 1.0/(k);
                    }
                }
                else if(gamma == Gamma.None) // none
                {
                    gammaArray = Dirichlet(gamma,k,N);
                    if(lambda_init == null)
                    {
                        lambda_init = new List<double>();
                        lambda_init.Add(H.Median(H.Divide(x,t)) / 2);
                        Console.WriteLine("Default gamma initialization, lamda: " + lambda_init.ToString());
                    }
                    for(int i = 0; i < A; i++)
                    {
                        mu[i]  = H.Multiply(lambda_init, (i+1))[0];
                        mix[i] = 1.0/(k);
                    }
                }
                else 
                {
                    gammaArray = Dirichlet(gamma,k,N);
                    if(lambda_init == null)
                    {
                        lambda_init = new List<double>();
                        lambda_init.Add(H.Median(H.Divide(x,t)));
                        Console.WriteLine("Default gamma initialization, lamda: " + lambda_init.ToString());
                    }
                    for(int i = 0; i < A; i++)
                    {
                        mu[i]  = H.Multiply(lambda_init, (i+1))[0];
                        mix[i] = 0.4/(k-1);
                    }
                    mix[0] = 0.6;
                }

                int[,] R = new int[(int)N, (int)A]; // should be populated with rand values btw 0 and 1
                int done = 0; 
                int retry2 = 0;
                int iter = 0; 
                double prev_likelihood = 1.0e20;
                int[] px;
                int[] pxi;


                while(done == 0 && retry2 == 0)
                {
                    iter = iter + 1;
                    done = 1;
                    if(addUniform == 1)
                    {
                        //px,pxi,outliers = Poisson_geometric(x,t,A,mix,mu,tries,N,alpha);
                        Dictionary<string, int[]> geo = Poisson_geometric(x,t,A,mix,mu,tries,N,alpha);
                        px = geo["px"];
                        pxi = geo["pxi"];
                        int outliers = geo["outliers"][0];
                    }
                    else 
                    {
                        Dictionary<string, int[]> poi = Poisson(x, t, A, mix, mu, tries, N);
                        px = poi["px"];
                        pxi = poi["pxi"];
                    }

                    if(retry2 == 0)
                    {
                        for(int i = 0; i < A; i++)
                        {
                            for(int contig = 0; contig < pxi[i,:].Count; contig++){
                                
                                if(px[contig] > 0)
                                {
                                    R[contig,i] = np.transpose((pxi[i,contig] * mix[i]) / px[contig]) 
                                }
                                else
                                {
                                    R[contig,i] = 0;
                                }
                            }                                
                            sumR[i]     = sum(R[:,i])

                            if(gamma != Gamma.None) // or gammaArray != null??
                            {
                                mix_new[i]  = ((1/N)*sumR[i] + (1/N)*(gammaArray[i]-1)) / (1 + (gammaArray.Sum-k)/N);
                            } 
                            else 
                            {
                                mix_new[i] = sumR[i] / N; 
                            }
                
                        }
                        
                        likelihood = sum(np.log(px+epsilon));

                        for(int j = 0; j < A; j++)
                        {
                            likelihood += sum(R[:,j]*np.log(mix[j]+epsilon));
                            if gamma is not None:
                              likelihood += (gamma[j]-1)*np.log(mix[j]+epsilon);
                        }

                        denomgeom                  = A*t*R[:,-1]*np.log(1-alpha/np.float64(t))
                        denomgeom[x-A*t*mu[0] < 0] = 0
                        denomgeom_sum              = sum(denomgeom)
                
                        denompois                     = 0


                            for i in range(0,int(A)-addUniform):
                    denompois += sum(R[:,i]*(i+1)*t)
                denom = denomgeom_sum + denompois
                
                numer = 0;
                for(int i = 0; i < (int)A-addUniform; i++)
                {
                    numer += sum(R[:,i]*x);
                }
                mu_new[0] = numer/denom;
                for(int i = 1; i = A; i++)
                {
                    mu_new[i] = mu_new[0]*(i+1);
                }

                alpha = sumR[-1] / sum(R[:,-1]*(x/np.float64(t)-mu[-1]+1));  // alpha geometric
                if(double.IsNaN(alpha))
                {
                    alpha = 10e-6;
                }
                    
                if(iter%100 == 0)
                {
                    //print '[%3d] L: %2.2f (change: %2.2e); sum (P(j|x)) = ; alpha = %2.2e; lambda = %2.2e' % 
                                 //(iter, likelihood, abs((likelihood - prev_likelihood)/likelihood),alpha,mu_new[0]),
                    for(int i = 0; i < A; i++)
                    {
                        //print '%2.2f ' % sumR[i],
                    }
                    
                    //print '; P(j) = ',
                    for(int i = 0; i < A; i++)
                    {
                        //print '%2.2f ' % mix[i],
                    }
                    
                    //print '\n',
                    //print 'Outliers:',outliers,
                    //print '\n'    
                }
                    
                done = (abs ((likelihood - prev_likelihood)/likelihood) < 1e-5); 
                
                if(iter >= max_iters)
                {
                    //print "Maximum number of iterations reached"
                    done = true;
                }

                if(done)
                {
                    //print "Number of iterations needed to converge: ", iter
                }
                
                prev_likelihood = likelihood;
                
                for(int i = 0; i < A; i++)
                {
                    mix[i]    = mix_new[i];
                    mu[i]     = mu_new[i];
                }
                
                    }
                }

            }
            return likelihood, mu[:k-addUniform], mix[:k-addUniform], alpha;
        }

        
        
       //Initialize Dirichlet parameter
       //INPUT:
       //gamma:    'haploid', 'diploid', 'uniform', or an list of gamma's of length n
       //n:        Number of models
       //I:        Hyperparameter that controls the prior's impact during model selection (I>0)
       
       //OUTPUT:
       //gamma array
         
        public float[] Dirichlet(Gamma gamma, int n, float I)
        {
            if(n < 2) { throw new Exception ("Poisson mixture needs at leat two Poisson models"); }
            if(I <= 0) { throw new Exception("Hyperparameter I should be larger then 0"); }

            float[] res = new float[n];
            for(int k = 0; k < n; k++)
            {
                res[k] = 1;
            }

            if(gamma == 1) // haploid
            {
                res[0] += I;
            } 
            else if(gamma == 2) // diploid
            {
                res[1] += I;
            }
            else if(gamma == 3) // allo-triploid
            {
                res[0] += (float)I/2;
                res[1] += (float)I/2;
            } 
            else if(gamma == 4) // di-triploid
            {
                res[1] += (float)I/2;
                res[2] += (float)I/2 ;       
            }
            //elif gamma is None or gamma == 'uniform':
            else if(gamma == -1 || gamma == 0) ///???
            {
                res[0] +=  (float)I/2;
                //res = (  np.ones(n)+(float)I/n  ); //???
            }
            else 
            {
                throw new Exception("Provide a valid gamma, see the docstring of em.Dirichlet");
            }

            return res;
        }

          private Dictionary<string, int[]> Poisson(int[] x, int[] t, float A, List<double> Mix, List<double> Mu, int tries, float N) 
        {
            //  # P(data|model i) and P(data)
            int[] pxi = new int[(int)N*(int)A];
            int[] px = new int[(int)N];

            for(int i = 0; i < A; i++)
            {
                  //# Fill in poisson
                  //pxi[i,:]          = stats.poisson.pmf(x,mu[i]*t)
                  //px                = px + mix[i] * pxi[i,:]
            }
            
            Dictionary<string, int[]> dict = new Dictionary<string, int[]>();
            dict.Add("px", px);
            dict.Add("pxi", pxi);
            return dict;
        }

           private Dictionary<string, int[]> Poisson_geometric(int[] x, int[] t, float A, List<double> mix, 
                List<double> mu, int tries, float N, double alpha)
        {
            int x_delay  = x - (int)(mu[mu.Count]* (int)t);
            double minxt    = mu[mu.Count];
            
            //# P(data|model i) and P(data)
            int[] pxi = new int[(int)N * (int)A];
            Array.Clear(pxi, 0, pxi.Length);

            int[] px = new int[N];
            Array.Clear(px, 0, px.Length);

            for(int i = 0; i < A; i++)
            {
                if(i < A-1)
                {
                    //# Fill in poisson
                    //pxi[i,:]          = stats.poisson.pmf(x,mu[i]*t) # Add machine precision
                    //px                = px + mix[i] * pxi[i,:]
                }
                else 
                {
                    //# Fill in geometric
                    //pxi[i,:] = stats.distributions.geom.pmf(x_delay,alpha/np.float64(t))
                    //px       = px + mix[i] * pxi[i,:]
                }
            }

            Dictionary<string, int[]> dict = new Dictionary<string, int[]>();
            dict.Add("px", px);
            dict.Add("pxi", pxi);
            dict.Add("outliers", np.sum([x_delay > 0]));
            return dict;
        }  
    

    }


        // helper class
        private static class H
        {
            public static decimal Median(List<int> xs)
            {
                xs.Sort();
                return xs[xs.Count / 2];
            }

            public static float Median(List<float> xs)
            {
                xs.Sort();
                return xs[xs.Count / 2];
            }

            // Replaces python division
            public static List<float> Divide(int[] a, int[] b)
            {
                List<float> c = new List<float>();
                for(int i = 0; i < a.Length; i++)
                {
                    c.Add((float)a[i]/(float)b[i]);
                }
                return c;
            }

            public static List<float> Divide(List<int> a, List<int> b)
            {
                List<float> c = new List<float>();
                for(int i = 0; i < a.Count; i++)
                {
                    c.Add((float)a[i]/(float)b[i]);
                }
                return c;
            }

            // multiply each of a by b
            public static List<double> Multiply(List<double> a, int b)
            {
                List<double> c = new List<double>();
                for(int i = 0; i < a.Count; i++)
                {
                    c.Add(a[i] * b);
                }
                return c;
            }
        }
    ////////////////////////////
    
        // Calculate CNV based on mixture model (PoissonMixture() below)
    private static class Cnv
    {
        private List<int> mixtures;
        private List<int> posteriors;
        private bool first;

        public Cnv(data, int nrsamples, int[] k, 
            double epsilon, max_iters, gammas, minSize, minCp, maxCP)
        {
            this.mixtures = new List<int>();
            this.posteriors = new List<int>();
            this.first = true;
        

        for(int i = 0; i < nrsamples; i++)
        {
            string figname = "poisson" + i + ".pdf";
            [counts,intervals] = GetCountsAndLengths(data, sample)
          [counts,intervals] = FilterContigs(counts,intervals,minSize,minCP,maxCP)

          Console.WriteLine("Fitting mixture model sample: "+i);
            this.mixtures.Add( new PoissonMixture(counts, intervals,k, epsilon,max_iters,gammas[sample]) );
        }

            
            //# Calculate optimal k base on a combined BIC
            //int[] BICall = zeros(len(self.mixtures[0].BICs));
            int[] BICall = new int[mixtures[0].Bics];
            Array.Clear(BICall, 0, BICall.Length);
        for sample in range(nrsamples):
            if self.mixtures[sample].BICs is not None:
              BICall += self.mixtures[sample].BICs
        best = argmin(BICall)
        print "Combined BIC optimal is: ", self.mixtures[0].ks[best]
        
        # Set new k
        for sample in range(nrsamples):
          self.mixtures[sample].setModelNumber(best)
        
        for sample in range(nrsamples):
            print "Calculating posteriors"
            self.posteriors.append(self.posteriorAllContigs(self.mixtures[sample],data,sample))

    }
        



          
    def printParameters(self,prefix):
      if prefix == "":
        modelfile = "model.txt"
      else:
        modelfile = prefix + ".model.txt"
      h = open(modelfile,"w")
      for mixture in self.mixtures:
        str = '%d\t%.6f\t%.6f' % (len(mixture.mu), mixture.mu[0], mixture.alpha)
        for i in range(len(mixture.mix)):
          str += '\t%.6f' % (mixture.mix[i])
        print >>h, str
      h.close()
      return
      

    def getCountsAndlenghts(self,data,sample):
      counts     = []; lengths    = []
      contigs = data.keys()
      for contig in contigs:
        counts.append(data[contig]['counts'][sample])
        lengths.append(data[contig]['clen'])
      return counts, lengths
    
    def cnv(self, sample1, sample2):
      pCnv = {}
      for node in self.posteriors[sample1]:
          pCnv[node] = self.probabilityCNV(self.posteriors[sample1][node],self.posteriors[sample2][node])
      return pCnv
  
    def cn(self,sample):
      res = {}
      for node in self.posteriors[sample]:
          res[node] = argmax(self.posteriors[sample][node]) + 1
      return res
        
    def isZeroCN(self, contigid, data, sample, multiplier):
      """Set CN to zero if it is lower than 'multiplier x lambda'"""
      mu     = self.mixtures[sample].mu[0]
      counts = data[contigid]['counts'][sample]
      if (counts/float64(data[contigid]['clen']) < (mu - multiplier*sqrt(mu))):
        return True
      else:
        return False
        
    #def printCN(self,am,sample,ref=None,file=stdout):
    #  if file != stdout:
    #      file = open(file,"w")
    #  cns = self.cn(sample)
    #  for node in am.ass.contigGraph.contigStats:
    #      count  = am.ass.contigGraph.contigStats[node]['counts'][0]
    #      length = am.ass.contigGraph.contigStats[node]['length']
    #      cn     = cns[node]
    #      if ref and ref.contigLocation.has_key(node):
    #          chr    = ref.contigLocation[node]['chr']
    #          start  = ref.contigLocation[node]['feat'].location.nofuzzy_start+1
    #          end    = ref.contigLocation[node]['feat'].location.nofuzzy_end
    #      else:
    #          chr    = ''
    #          start  = 0
    #          end    = 0
    #      # contig ID, length, count, CN, CHR, start, stop
    #      print >>file,'%s\t%d\t%d\t%d\t%s\t%d\t%d' % (node,length,count,cn+1,chr,start,end)
    #  if file != stdout:
    #      file.close()    
    
    def posteriorAllContigs(self,mixture, data, sample):
      """Calculate posteriar probabilities for each contig"""
      prob = {}
      for node in data.keys():
          # TODO fix this in absence of am
          #if am.ass.contigGraph.contigStats[node]['sampleMems'][sample]:
          if data[node]['counts'][sample] >= 1:
            prob[node] = mixture._posterior(data[node]['counts'][sample],data[node]['clen'])
          else:
            prob[node] = [0]*len(mixture.mix)
      return prob
    
    def probabilityCNVAllContigs(self, s1, s2):
      res   = {}
      post1 = self.posteriors[s1]
      post2 = self.posteriors[s2]
      for node in post1.keys():
          res[node] = self.probabilityCNV(post1[node], post2[node])
      return res
        
    def probabilityCNV(self, post1, post2):
      """Probability a contig has a different copy number in the two samples
      Input: Two lists with posterior probabilities for all copy numbers (k)"""
      p_equal = 0
      for k in range(len(post1)):
          p_equal += post1[k] * post2[k]
      return 1-p_equal
        
    def _filterContigs(self, counts, intervals, minSize, minCP, maxCP):
      """Filter contigs with too little reads or too small
      minSize: minimum contig size
      minCP  : minimum number of reads per 100 nucleotides [1]
      maxCP  : maximum number of reads per 100 nucleotides [mean(count/100)*5]"""
      counts    = array(counts) 
      intervals = array(intervals)
      countsCP  = (counts / float64(intervals)) * 100
      if maxCP == "Default":
          maxCP = mean(countsCP) * 5
          
      print "Contigs < ", minSize, " are not used to fit the model"
      
      ind  =  (intervals >= minSize) & (countsCP  >= minCP) & (countsCP <= maxCP)                 
      return counts[ind], intervals[ind]





    
    
}












              // for given cluster, run
        private void RunMagnolya()
        {
            //magnolya -i counts.txt -p ecoli -g haploid \
            //-l contiglocs.txt -m 5 -s 500 -S 500
        }


        private static class CnvCall
        {

        }
*/