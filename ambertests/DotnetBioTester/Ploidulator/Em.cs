/*using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Bio.Algorithms.Metric
{
    //Expectation Optimization procedure used by Magnolya
    public class Em
    {

        private void PoissonGeometric(int x, int t, int A, List<int> mix, List<int> mu, int tries, int N, int alpha) 
         
        {
            int x_delay  = x - (mu[mu.Count - 1]*t);
            int minxt    = mu[mu.Count - 1];

            //# P(data|model i) and P(data)
            //pxi = np.array([np.zeros(N)]*int(A)); 
            //px = np.zeros(N);
            int[] pxi = new int[N * A];
            int[] px = new int[N];

            for(int i = 0; i < A; i++)
            {
                if(i < A-1)
                {
                    //Fill in poisson
                    //pxi[i,:]          = stats.poisson.pmf(x,mu[i]*t) # Add machine precision
                    //px                = px + mix[i] * pxi[i,:]
                    px = px + mix[i] * pxi[i,:];
                } 
                else 
                {
                    //Fill in geometric
                    //pxi[i,:] = stats.distributions.geom.pmf(x_delay,alpha/np.float64(t))
                    //px       = px + mix[i] * pxi[i,:]
                }
            }
            //return px, pxi, np.sum([x_delay > 0])    
            return new Tuple<int[],int[], int>(px, pxi, sum([x_delay > 0]));
        }


        private Dictionary<string, int[]> Poisson(int x, int t, int A, List<int> mix, List<int> mu, int tries, int N) 
        {
            //  # P(data|model i) and P(data)
            int[] pxi = new int[N*A];
            int[] px = new int[N];

            for(int i = 0; i < A; i++)
            {
                  //# Fill in poisson
                  pxi[i,:]          = stats.poisson.pmf(x,mu[i]*t)
                  px                = px + mix[i] * pxi[i,:]
            }
            
            Dictionary<string, int[]> dict = new Dictionary<string, int[]>();
            dict.Add("px", px);
            dict.Add("pxi", pxi);
            return dict;
        }


        
       //Initialize Dirichlet parameter
       //INPUT:
       //gamma:    'haploid', 'diploid', 'uniform', or an list of gamma's of length n
       //n:        Number of models
       //I:        Hyperparameter that controls the prior's impact during model selection (I>0)
       
       //OUTPUT:
       //gamma array
         
        public float[] Dirichlet(int gamma, int n, float I)
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

   

        // k = range of cluster numbers in a list (renamed cluster)
        public Tuple<List<List<int>>, List<double>> EmPoissonBic(List<int> count, List<int> interval, 
                List<int> cluster = null, double epsilon = 2e-13, 
                int max_iters = 100, int gamma = 2, string figname = "poisson.pdf", int addUniform = 1,
                List<double> initPercLists = null, bool bw = true, int plotContigSize = 1000, int maxVal = int.MaxValue, List<int> ymax = null)
        {
            if(count.Count == 0) { Console.WriteLine("Warning: Zero reads detected"); return null; }

            Console.WriteLine("****Starting model selection using BIC****");
            Console.WriteLine("Number of contigs used to train EM: "+count.Count);
            Console.WriteLine("Median contig length: ", Median(interval));
    
            cluster = cluster == null ? new List<int>(new int[]{3,4,6,8,10,12,15,20}) : cluster;
            
            List<int> likelihoods = new List<int>();    List<int> mus = new List<int>();
            List<int> mixs = new List<int>();           List<int> alphas = new List<int>();

            // for each cluster
            foreach(int sequence in cluster)
            {
                Dictionary<string, int> results = EmPoissonInit(count, interval, 
                        sequence, 0, epsilon, max_iters, gamma, 
                        figname, addUniform, initPercLists, plotContigSize, maxVal);
                likelihoods.Add(results["likelihood"]);
                mus.Add(results["mu"]);
                mixs.Add(results["mix"]);
                alphas.Add(results["alpha"]);
            }

            //# Calculate Bayesian Information Criterion
            int n = count.Count;
            List<double> BICS = new List<double>();

            for(int i = 0; i < likelihoods.Count; i++)
            {
                double val = (-2 * likelihoods[i] + (cluster[i]+1+2*addUniform) * Math.Log(n)); // # pi's + lambda + geom (alpha and pi);
                BICS.Add(val);
                Console.WriteLine(String.Format("k = %d; BIC = %2.2f; L = %2.2f; lambda = %2.2e;",
                          cluster[i], BICS[i], likelihoods[i], mus[i])); // was mus[i][0]    
            }
            int best = BICS.IndexOf(BICS.Min());
            
            List<List<int>> retList = new List<List<int>>(new List<int>[]{likelihoods, mus, mixs, cluster, alphas});
            return new Tuple<List<List<int>>, List<double>>(retList, BICS);
        }
    

        public Dictionary<string, int> EmPoissonInit(List<int> count, List<int> interval, 
                int k = 10, int plot_mixture = 0, double epsilon = 2e-13, int max_iters = 100, 
                int gamma = 1, string figname = "poisson.pdf", int addUniform = 1, 
                List<double> initPercList = null, int plotContigSize = 1000, int maxval = int.MaxValue)
        {
            //INITPERCLIST    List with percentiles to initialize the EM
            if(initPercList == null && (gamma != 1 && gamma != 2))
            {
                initPercList = new List<double>(new double[]{1, 5, 12.5, 25, 50, 70});   
                Console.WriteLine("Performing initializations at 1, 5, 12.5, 25 and 50th percentile");
            } 
            else if(initPercList == null && (gamma == 1 || gamma == 2))
            {
                initPercList = new List<double>(new double[]{50});   
                Console.WriteLine("Performing initialization at 50th percentile");
            }
            
            List<int> likelihoods = new List<int>();    List<int> mus = new List<int>();
            List<int> mixs = new List<int>();           List<int> alphas = new List<int>();

            foreach(double perc in initPercList) 
            {
                // lambda for 1, 5, 12.5, 25, 50, 70 is calculated below (currently returns null though)
                List<double> lambda_init = ScoreAtPercentile(Divide(count, interval), perc); // todo fixme

                List<int> results = EmPoisson(count,interval,k,0,epsilon,max_iters,
                        gamma, figname, addUniform, lambda_init, plotContigSize, maxval);
                int likelihood = results[0];
                int mu = results[1];
                int mix = results[2];
                int alpha = results[3];

                likelihoods.Add(likelihood);
                mus.Add(mu);
                mixs.Add(mix);
                alphas.Add(alpha);
                
                Console.WriteLine(String.Format("Lambda init: %2.2e for %2.2f-th percentile -> L = %2.2f  ",
                          lambda_init, perc,likelihood)); 
            }
            int best = likelihoods.IndexOf(likelihoods.Max());
            
            Dictionary<string, int> returnDict = new Dictionary<string, int>();   
            returnDict.Add("likelihood", likelihoods[best]);
            returnDict.Add("mu", mus[best]);
            returnDict.Add("mix", mixs[best]);
            returnDict.Add("alpha", alphas[best]);
            return returnDict;
        }
    
    
             // todo fixme n ote that there is a much faster alternative
        decimal Median(List<int> xs)
        {
            xs.Sort();
            return xs[xs.Count / 2];
        }

        float Median(List<float> xs)
        {
            xs.Sort();
            return xs[xs.Count / 2];
        }

        // Replaces python division
        private List<float> Divide(int[] a, int[] b)
        {
            List<float> c = new List<float>();
            for(int i = 0; i < a.Length; i++)
            {
                c.Add((float)a[i]/(float)b[i]);
            }
            return c;
        }

        private List<float> Divide(List<int> a, List<int> b)
        {
            List<float> c = new List<float>();
            for(int i = 0; i < a.Count; i++)
            {
                c.Add((float)a[i]/(float)b[i]);
            }
            return c;
        }

        // multiply each of a by b
        private List<double> Multiply(List<double> a, int b)
        {
            List<double> c = new List<double>();
            for(int i = 0; i < a.Count; i++)
            {
                c.Add(a[i] * b);
            }
            return c;
        }

        private List<double> ScoreAtPercentile(List<float> data, double per)
        {
            if(per < 0 || per > 100) { throw new Exception("The percentile should be between 0. and 100. !"); }
            
            /*private void _quantiles1D()
            {
                        x = np.sort(data.compressed())
                n = len(x)
                if n == 0:
                    return ma.array(np.empty(len(p), dtype=float), mask=True)
                elif n == 1:
                    return ma.array(np.resize(x, p.shape), mask=nomask)
                aleph = (n*p + m)
                k = np.floor(aleph.clip(1, n-1)).astype(int)
                gamma = (aleph-k).clip(0,1)
                return (1.-gamma)*x[(k-1).tolist()] + gamma*x[k.tolist()]
            }
             

    //# Initialization & checks ---------
    data = ma.array(a, copy=False)
    if data.ndim > 2:
        raise TypeError("Array should be 2D at most !")
    #
    if limit:
        condition = (limit[0]<data) & (data<limit[1])
        data[~condition.filled(True)] = masked
    #
    p = np.array(prob, copy=False, ndmin=1)
    m = alphap + p*(1.-alphap-betap)
    # Computes quantiles along axis (or globally)
    if (axis is None):
        return _quantiles1D(data, m, p)
    return ma.apply_along_axis(_quantiles1D, axis, data, m, p)
            

            return null;
        }

    
      //   * [L, MU, S2, MIX] = EM (COUNT, INTERVAL, K, STEPK, FIX_S2, REG, PLOT, EPS, MAXITER)
    
    //Uses the EM algorithm to fit a mixture of K 1D Poissons, where the means
    //and the variances of the models are spaced as (1 + k/STEPK), k = 0..K-1.
    //Optional parameters are: 
    //FIX_S2	set the variance to the mean (1, dflt.) or estimate (0)
    //REG     regularisation constant (dflt. 0)
    //PLOT    plot progress (dflt. 0)
    //EPS     epsilon, the change in likelihood to stop training (dflt. 1e-5)
    //MAXITER maximum number of iterations (dlft. 100)
    //Returns the likelihood L, the means MU(), the variances S2() and priors MIX.
         
        public List<int> EmPoisson(List<int> count, List<int> interval, int k = 10, int plot_mixture = 0, 
            double epsilon = 2e-13, int max_iters = 100, int gamma = 1, string figname = "poisson.pdf", 
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
                        
            double min_data = (Divide(x,t).Min()) > epsilon ? Divide(x,t).Min() : epsilon;
            float max_data = Divide(x,t).Max(); //# Normalised for interval
            int tries = 0; 
            int retry = 1;
            float[] gammaArray;


            while(retry == 1 && tries <= 5)
            {
                // todo fixme i don't think these are all right
                if(gamma == 1) // haploid
                {
                    gammaArray = Dirichlet(gamma,k,N);
                    if(lambda_init == null)
                    {
                        lambda_init = new List<double>(); // should just be double
                        lambda_init.Add(Median(Divide(x,t)));
                        Console.WriteLine(String.Format("Haploid initialization, lamda: " , lambda_init.ToString()));
                    }
                    for(int i = 0; i < A; i++)
                    {
                        mu[i]  = Multiply(lambda_init, (i+1))[0];
                        mix[i] = 0.4/(k-1);
                    }            
                    mix[0] = 0.6;
                }
                else if(gamma == 2) // diploid
                {
                    gammaArray = Dirichlet(gamma,k,N);
                    if(lambda_init == null)
                    {
                        lambda_init = new List<double>();
                        lambda_init.Add(Median(Divide(x,t)) / 2);
                        Console.WriteLine(String.Format("Diploid initialization, lamda: " , lambda_init.ToString()));
                    }
                    for(int i = 0; i < A; i++)
                    {
                        mu[i]  = Multiply(lambda_init, (i+1))[0];
                        mix[i] = 0.4/(k-1); // check whether this shoul be same as above or diff
                    }            
                    mix[0] = 0.6;
                }
                else if(gamma == 3) // allo-triploid
                {
                    gammaArray = Dirichlet(gamma,k,N);
                    if(lambda_init == null)
                    {
                        lambda_init = new List<double>();
                        lambda_init.Add(Median(Divide(x,t)) / 3);
                    }
                    for(int i = 0; i < A; i++)
                    {
                        mu[i]  = Multiply(lambda_init, (i+1))[0];
                        mix[i] = 0.2/(k-2);
                    }
                    mix[1] = 0.4;
                    mix[2] = 0.4;   
                    Console.WriteLine("haplo-diploid initialization, lamda: " + lambda_init.ToString());
                }
                else if(gamma == 4) // di-triploid
                {
                    gammaArray = Dirichlet(gamma,k,N);
                    if(lambda_init == null)
                    {
                        lambda_init = new List<double>();
                        lambda_init.Add(Median(Divide(x,t)) / 3);
                        Console.WriteLine("Di-triploid initialization, lamda: " + lambda_init.ToString());
                    }
                    for(int i = 0; i < A; i++)
                    {
                        mu[i]  = Multiply(lambda_init, (i+1))[0];
                        mix[i] = 0.2/(k-2);
                    }
                    mix[1] = 0.4;
                    mix[2] = 0.4;   
                }
                else if(gamma == 5) // uniform
                {
                    gammaArray = Dirichlet(gamma,k,N);
                    if(lambda_init == null)
                    {
                        lambda_init = new List<double>();
                        lambda_init.Add(Median(Divide(x,t)) / 3);
                        Console.WriteLine("Default initialization, lamda: " + lambda_init.ToString());
                    }
                    for(int i = 0; i < A; i++)
                    {
                        mu[i]  = Multiply(lambda_init, (i+1))[0];
                        mix[i] = 1.0/(k);
                    }
                }
                else if(gamma == -1) // none
                {
                    gammaArray = Dirichlet(gamma,k,N);
                    if(lambda_init == null)
                    {
                        lambda_init = new List<double>();
                        lambda_init.Add(Median(Divide(x,t)) / 2);
                        Console.WriteLine("Default gamma initialization, lamda: " + lambda_init.ToString());
                    }
                    for(int i = 0; i < A; i++)
                    {
                        mu[i]  = Multiply(lambda_init, (i+1))[0];
                        mix[i] = 1.0/(k);
                    }
                }
                else 
                {
                    gammaArray = Dirichlet(gamma,k,N);
                    if(lambda_init == null)
                    {
                        lambda_init = new List<double>();
                        lambda_init.Add(Median(Divide(x,t)) / 2);
                        Console.WriteLine("Default gamma initialization, lamda: " + lambda_init.ToString());
                    }
                    for(int i = 0; i < A; i++)
                    {
                        mu[i]  = Multiply(lambda_init, (i+1))[0];
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
                        Dictionary<string, int[]> geo = Poisson_geometric(x,t,A,mix,mu,tries,N,alpha);
                        px = geo["px"];
                        pxi = geo["pxi"];
                        int outliers = geo["outliers"][0];
                    }
                    else 
                    {
                        Dictionary<string, int[]> poi = Poisson(x,t,A,mix,mu,tries,N);
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
                            if gamma is not None:
                                mix_new[i]  = ((1/N)*sumR[i] + (1/N)*(gamma[i]-1)) / (1 + (sum(gamma)-k)/N);
                            else:
                                mix_new[i] = sumR[i] / N;  # Gemiddelde aandeel alle datapunten 
                
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

                //alpha = sumR[sumR.Count] / sum(R[:,-1]*(x/np.float64(t)-mu[-1]+1));  // alpha geometric
                alpha = sumR[sumR.Count-1] / sum(  R[:,-1]  *  (H.Divide(x,t)  - mu[mu.Count-1]+1)  );  // alpha geometric
                if(double.IsNaN(alpha))
                {
                    alpha = 10e-6;
                }
                    
                if(iter%100 == 0)
                {
                    Console.Write(String.Format("[%3d] L: %2.2f (change: %2.2e); sum (P(j|x)) = ; alpha = %2.2e; lambda = %2.2e"), 
                        iter, likelihood, abs((likelihood - prev_likelihood)/likelihood),alpha,mu_new[0]);
                    for(int i = 0; i < A; i++)
                    {
                        Console.Write(sumR[i]); //print '%2.2f ' % sumR[i],
                    }
                    Console.Write("; P(j) = ");
                    for(int i = 0; i < A; i++)
                    {
                        Console.Write(mix[i]); //print '%2.2f ' % mix[i],
                    }
                    Console.Write("\nOutliers: "+outliers + "\n");
                }
                    
                done = (abs ((likelihood - prev_likelihood)/likelihood) < 1e-5); 
                
                if(iter >= max_iters)
                {
                    Console.WriteLine("Maximum number of iterations reached");
                    done = true;
                }

                if(done)
                {
                    Console.WriteLine("Number of iterations needed to converge: "+iter);
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
}
*/