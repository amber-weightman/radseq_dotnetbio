using Accord.Statistics.Distributions.Univariate;
using Accord.Statistics.Models.Markov;
using Accord.Statistics.Models.Regression;
using Accord.Statistics.Models.Regression.Fitting;
using Bio;
using Bio.Algorithms.Alignment;
using Bio.Algorithms.Metric;
using Bio.IO.BAM;
using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Media;
using System.Windows.Shapes;
using System.Windows.Threading;

namespace Ploidulator
{
    /// <summary>
    /// An IMetric calculates various metric values from a list of SAMAlignedSequences, where each
    /// list of sequences belongs to a cluster. Clusters are determined based on sequence similarity, 
    /// e.g. as calculated by an MCL graph clustering algorithm, or any other clustering approach
    /// </summary>
    public class ClusterMetricPloidulator : IMetric
    {

        #region Private Fields

        private int expectedPloidy;

        /// <summary>
        /// The number of samples represented in the entire population
        /// </summary>
        private int numSamples;
        
        private List<Tuple<char, char>> genotype = null;


        // buckets of everything
        private List<SAMAlignedSequence> readsInPloidyForIndividuals = null;
        private List<SAMAlignedSequence> readsNotInPloidyForIndividuals = null;
        int totQualInPloidy = 0;
        int totQualNotInPloidy = 0;

        // 
        private Dictionary<string, List<SAMAlignedSequence>> readsInPloidyForIndividualsDict = null;
        private Dictionary<string, List<SAMAlignedSequence>> readsNotInPloidyForIndividualsDict = null;

        /// <summary>
        /// A dictionary where each distinct query sequence is represented as a key, and 
        /// a list of all individual sequences that share the same query string is stored
        /// as the value
        /// </summary>
        private Dictionary<String, List<SAMAlignedSequence>> sequenceDict = null;

        /// <summary>
        /// A dictionary where each sample individual is represented as a key, and 
        /// a list of all individual sequences for that individual is stored
        /// as the value
        /// </summary>
        private Dictionary<String, List<SAMAlignedSequence>> sampleDict = null;

        /// <summary>
        /// Experimental
        /// </summary>
        private Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> sampleSequenceDict = null;

        /// <summary>
        /// Experimental
        /// </summary>
        private Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> sequenceSampleDict = null;

        /// <summary>
        /// List of all sequences.
        /// </summary>
        private List<SAMAlignedSequence> sequences = null;

        //private SequenceAlignmentMap sequenceMap = null;

        #endregion

        #region Constructors

        /// <summary>
        /// The default constructor.
        /// </summary>
        public ClusterMetricPloidulator()
        {
            //throw new NotImplementedException();
        }

        public ClusterMetricPloidulator(int ploidy, int samples)
        {
            expectedPloidy = ploidy;
            numSamples = samples;
            
        }

        #endregion

        #region Properties

        /// <summary>
        /// Cluster ID (unique ID of reference sequence against which all sequences in the cluster
        /// are aligned)
        /// (If sequences currently being handled are from multiple clusters, only the first Id
        /// will be returned)
        /// </summary>
        public string Id { get { return id; } }
        private string id;

        /// <summary>
        /// Total number of all sequences in the cluster
        /// </summary>
        public int CountAll { get { return countAll; } }
        private int countAll;

        /// <summary>
        /// Total number of distinct sequences in the cluster
        /// </summary>
        public int CountDistinct { get { return countDistinct; } }
        private int countDistinct;

        /// <summary>
        /// Total number of samples (individuals) represented in the cluster
        /// </summary>
        public int CountSamples { get { return countSamples; } }
        private int countSamples;

        /// <summary>
        /// A dictionary where each distinct query sequence is represented as a key, and 
        /// a list of all individual sequences that share the same query string is stored
        /// as the value
        /// </summary>
        public Dictionary<String, List<SAMAlignedSequence>> SequenceDict { get { return sequenceDict; } }

        /// <summary>
        /// A dictionary where each sample individual is represented as a key, and 
        /// a list of all individual sequences for that individual is stored
        /// as the value
        /// </summary>
        public Dictionary<String, List<SAMAlignedSequence>> SampleDict { get { return sampleDict; } }

        /// <summary>
        /// List of all sequences.
        /// </summary>
        public List<SAMAlignedSequence> Sequences { get { return sequences; } }

        /// <summary>
        /// Sequences as map
        /// </summary>
        //public SequenceAlignmentMap SequenceMap { get { return sequenceMap; } }

        /// <summary>
        /// Experimental.
        /// </summary>
        public List<int> FrequencyDistributionSequences { get { return frequencyDistributionSequences; } }
        private List<int> frequencyDistributionSequences;

        /// <summary>
        /// Experimental.
        /// </summary>
        public List<int> FrequencyDistributionSamples { get { return frequencyDistributionSamples; } }
        private List<int> frequencyDistributionSamples;

        /// <summary>
        /// The total number of sequences which are within ploidy (i.e. the top two represented query
        /// sequences)
        /// </summary>
        //public int CountInPloidy { get { return GetCountInPloidy(2, GetFrequencyDistribution(SequenceDict)); } }

        /// <summary>
        /// The total number of sequences which are beyond ploidy (i.e. total count excluding the 
        /// top two represented query sequences)
        /// </summary>
        //public Double BeyondPloidy { get { return GetPercentBeyondPloidy(2, GetFrequencyDistribution(SequenceDict)); } }

        /// <summary>
        /// Average frequencies for each sequence, calculated per sample per cluster and averaged
        /// out to cluster
        /// </summary>
        public double[] ClustSeqFrequencies { get { return clustSeqFrequencies; } }

        /// <summary>
        /// Some measure of alignment quality using MAPQ
        /// </summary>
        public double[] ClustAlignmentQualities { get { return clustAlignmentQualities; } }

        /// <summary>
        /// Placeholder (determine if cluster is good based on metrics etc)
        /// if the metric determines that we want to keep this read for downstream analysis
        /// </summary>
        public bool Good { get { return good; } set { this.good = value; } }

        
        private bool good = false;

        /// <summary>
        /// Placeholder (determine if cluster is bad based on metrics etc)
        /// </summary>
        public bool Bad { get { return false; } }


        public double Dirt { get { return dirt; } }
        private double dirt;

        /// <summary>
        /// being the avg qual for distinct reads per indiv which are within ploidy (the reads on the right side of cluster dirt)
        /// </summary>
        public double AlignmentQual { get { return alignmentQual; } }
        private double alignmentQual;

        /// <summary>
        /// Are all individuals represented in the cluster? This is the percentage of whole-population representation
        /// </summary>
        public double PopulationPercentage { get { return populationPercentage; } }
        private double populationPercentage;

        public double[] SampleReadCountsAll { get { return sampleReadCountsAll; } }
        private double[] sampleReadCountsAll;

        public double[] SampleReadCountsDistinct { get { return sampleReadCountsDistinct; } }
        private double[] sampleReadCountsDistinct;

        public double[] AlignmentQualities { get { return alignmentQualities; } }
        private double[] alignmentQualities;

        public double[] ReadQualities { get { return readQualities; } }
        private double[] readQualities;

        #endregion

        #region Public Methods

        /// <summary>
        /// String of tab-separated values for writing to file by MetricFormatter.
        /// todo fixme we want the MetricFormater to decide how to write each line, but
        /// for now until I know what values are being written I will let the Metric handle this.
        /// </summary>
        public string ToFileString() // print file string
        {
            // PRINT THIS OUT TO SHOW DAD
            // for each in samplereadcountsAll, i want to know how much it deviates from the average
            // this requires some sort of frequency distribution
            // the array of nums that applies for fx dist is below
            /*string ss = "\nSamplereadcountsAll: ";
            foreach (double d in SampleReadCountsAll) { ss += (d + "\t"); }
            ss += ("\navg: " + SampleReadCountsAll.Average());*/

            // MAYBE ALSO SHOW DAD PRINTOUT FOR clustSeqFrequencies

            // todo update
            string header = (this.Id == "0") ? "cluster_id\tcount_all_reads\tcount_distinct_reads\tnum_individuals\tdirt\talignment_qualities_all:in:out\tploidy_disagreement_unnormalised\tread_qualities_all:in:out\tpopulation_percentage" : "";

            return header + Environment.NewLine + Id + "\t" + CountAll + "\t" +
                CountDistinct + "\t" +
                CountSamples + 
                "\t" + Dirt + "\t" + AlignmentQualities[0] + " : " + AlignmentQualities[1] + " : " + AlignmentQualities[2] + "\t(" + ploidyDisagreement + ")\t" +
                ReadQualities[0] + " : " + ReadQualities[1] + " : " + ReadQualities[2] + "\t" + PopulationPercentage + "\t" + Math.Round(SampleReadCountsAll.Average(), 2)
                //+ ss + Environment.NewLine
                ;
            ;
        }
        /// <summary>
        /// List of frequencies for each cluster, from most to least frequent (sums to 1). Ideally top [expectedPloidy] sequences
        /// would together be ~1 (100%)
        /// </summary>
        private double[] clustSeqFrequencies = null;

        /// <summary>
        /// The average alignment quality of each (distinct??) read in this cluster
        /// </summary>
        private double[] clustAlignmentQualities = null;



        /// <summary>
        /// Calculate metric values from the given list of sequences.
        /// </summary>
        /// <param name="sequences">Sequences to add.</param>

        public void Calculate(List<SAMAlignedSequence> sequences)
        {
            Console.Write(Id+"-");
            // Create various structures to store or index the sequence data in different ways
            this.sequences = sequences;
            sequenceDict = MakeSequenceDict(sequences);
            sampleDict = MakeSampleDict(sequences);
            sampleSequenceDict = MakeNestedSequenceDict(sampleDict);
            sequenceSampleDict = MakeNestedSampleDict(sequenceDict);
            
            // Get the counts once and once only
            countAll = sequences.Count;
            countDistinct = sequenceDict.Count;
            countSamples = sampleDict.Count;
            
            
            // Iterate through each structure the minimum number of times to calculate as many values as we can from it
            Console.Write("a");
            IterateSequences();
            Console.Write("b");
            IterateSequenceDict();
            Console.Write("c");
            IterateSampleDict();
            Console.Write("d");
            IterateSequenceSampleDict();
            Console.Write("e");
            IterateSampleSequenceDict();
            Console.Write("f");

            // Simple things that do not need to iterate to get their values
            populationPercentage = Math.Round(CountSamples / (double)numSamples, 2);
            Console.Write("g");

            // At this point every value should be set and nothing more should need to be done
            readQualities = FindReadQualities(); // takes a long time
            Console.Write("h");

            // CHECK WHICH DIRT CALC IS MORE ACCURATE
            dirt = Math.Round(1 - GetCountInPloidy(expectedPloidy, clustSeqFrequencies), 2);
            //double newDirt = readsNotInPloidyForIndividualsDict.Count / (double)(readsInPloidyForIndividualsDict.Count + readsNotInPloidyForIndividualsDict.Count);



            Console.Write("i");
            

            //clustAlignmentQualities = FindAlignmentQualities(sequenceDict); // quality; quantity

            Console.WriteLine(this.ToFileString()); // todo fixme this line of course should be removed
        }

        #region iterateSequences
        private void IterateSequences()
        {
            // set id (at the moment this doesn't iterate but in the future it should)
            id = SetId(sequences); 
                //double rqAll = FindReadQualities(sequences); (OCCURS BUT I HAVE IGNORED IT FOR NOW)
                // base frequencies (UNUSED BUT HERE)
            
            foreach (List<SAMAlignedSequence> seqList in sampleDict.Values)
            {
                Dictionary<char, double>[] fxx = BaseFrequencies(seqList, (double)expectedPloidy);
            }

            Dictionary<char, double>[] fx = BaseFrequencies(sequences, (double)expectedPloidy + 1);



            /*Tuple<char, char>[] geno = new Tuple<char, char>[GetSequence(sequences[0]).Length];

            for (int i = 0; i < geno.Length; i++ )
            {
               // AddGenotype(geno[i], seq.ch);
                
            }*/
            
        }

        private class BasePair<T1, T2>
        {
            public T1 First { get; set; }
            public T2 Second { get; set; }
        }

      

        private string SetId(List<SAMAlignedSequence> sequences)
        {
            if (!sequences[0].Flag.HasFlag(SAMFlags.UnmappedQuery))
            {
                String rname = sequences[0].RName;
                return (sequences != null) ? rname : null;
            }
            else
            {
                throw new Exception("Unmapped query exception"); // todo particular kind of exception?
            }
        }

        // THIS METHOD IS CURRENTLY UNUSED BUT THERE MAY BE A PURPOSE FOR IT LATER
        private Dictionary<char, double>[] BaseFrequencies(List<SAMAlignedSequence> seqs, double dirtCutoff)
        {
            double totalGDirt = 0;
            // length of REFERENCE
            Dictionary<char, double>[] freqList = new Dictionary<char, double>[GetSequence(seqs[0]).Length];

            //For each position, get each nucleotide and its relative frequency
            foreach (SAMAlignedSequence seq in seqs)
            {
                string seqStr = GetSequence(seq);
                int count = 0;
                foreach (char c in seqStr.ToUpper().ToCharArray())
                {
                    if (freqList[count] != null && freqList[count].ContainsKey(c))
                    {
                        ++freqList[count][c];
                    }
                    else
                    {
                        if (freqList[count] == null)
                        {
                            freqList[count] = new Dictionary<char, double>();
                        }
                        freqList[count].Add(c, 1);
                    }
                    count++;
                }
            }

            // for each position
            for (int i = 0; i < freqList.Length; i++ )
            {
                double numBases = freqList[i].Values.Sum();
                char[] chars = freqList[i].Keys.ToArray();
                foreach (char c in chars)
                {
                    freqList[i][c] = Math.Round((freqList[i][c] / numBases), 2);
                    if (freqList[i][c] == 0) { freqList[i].Remove(c); }
                }

                var sorted = (from b in freqList[i]
                                          orderby b.Value descending
                                          select b)
                    .ToDictionary(pair => pair.Key, pair => pair.Value);
                freqList[i] = sorted;



                int j = 0;
                double gDirt = 0;

                foreach (char c in freqList[i].Keys)
                {
                    //Console.Write("(" + c + ", " + freqList[i][c] + "), ");
                    if(j++ >= dirtCutoff)
                    {
                        gDirt += freqList[i][c];
                        
                    }

                }
                if(gDirt > 0)
                {
                    totalGDirt += gDirt;
                    //Console.Write("\t[" + gDirt + "] ");
                }
                //Console.WriteLine("");
                

            }
            //Console.WriteLine("--------------------" + totalGDirt + "--------------------");
            
            
            return freqList;
        }

        #endregion

        #region anySequenceList


        #endregion




        #region IterateSequenceDict

        // todo fixme this fails on clust 101
        private void IterateSequenceDict()
        {
            /*int i = 0;
            foreach(KeyValuePair<string, List<SAMAlignedSequence>> kvp in sequenceDict){
            }*/
            Console.Write("A");
            frequencyDistributionSequences = GetFrequencyDistribution(SequenceDict);
            Console.Write("B");
            alignmentQualities = FindAlignmentQualities(/*SequenceDict*/); // quality; quantity
            Console.Write("C");
            alignmentQual = alignmentQualities[1];
            Console.Write("D");
        }

        // this is used for either seq or sample dict
        private List<int> GetFrequencyDistribution(Dictionary<String, List<SAMAlignedSequence>> dict)
        {
            // todo aw probably a more efficient way
            List<int> frequencies = new List<int>();
            foreach (KeyValuePair<String, List<SAMAlignedSequence>> seq in dict)
            {
                frequencies.Add(seq.Value.Count);
            }
            return frequencies;
        }

        /// <summary>
        /// alignment quality rating PER DISTINCT SEQUENCE (only used on seq dict)
        /// dict sorted by sequence frequency within cluster
        /// i have assumed for now that this is only ever used with sequenceDict
        /// ploidy is per individual
        /// // this is skewed where a sequence appears in both in and out of ploidy buckets
        /// 
        /// 255 indicates qual is not available
        /// 29 is max
        /// read qual 0 is bad.
        /// </summary>
        /// <param name="dict"></param>
        /// <returns></returns>
        private double[] FindAlignmentQualities()
        {
            int i;
            int[] quals = new int[CountDistinct]; // one for every sequence in the map, just get its qualities
            double[] qualsIn = new double[readsInPloidyForIndividualsDict.Count]; // one for every seq which is in ploidy
            double[] qualsOut = new double[readsNotInPloidyForIndividualsDict.Count]; // one for every sequence which is out of ploidy

            Debug.Assert(CountAll == (readsNotInPloidyForIndividuals.Count + readsInPloidyForIndividuals.Count));
            Debug.Assert(CountDistinct <= (qualsIn.Length + qualsOut.Length));
            // second part of this will have more if some individuals have diff top 2 sequences from other individuals
            // disagreement is the number of distinct reads that are in the out of ploidy bucket despite being top-ploidy reads for some individusals
            // could be indicative of those individuals being 'triploid'?
            // the number should be small
            ploidyDisagreement = Math.Round(((qualsIn.Length + qualsOut.Length) - CountDistinct) / (double)CountDistinct, 2);
            //Console.WriteLine("ploidy disagreement " + ploidyDisagreement);

            i = 0;
            foreach (List<SAMAlignedSequence> seqList in sequenceDict.Values)
            {
                quals[i++] = seqList[0].MapQ; // assuming all identical sequences have the same quality reading
            }

            i = 0;
            foreach (List<SAMAlignedSequence> seqList in readsInPloidyForIndividualsDict.Values)
            {
                qualsIn[i++] = seqList[0].MapQ;
            }

            i = 0;
            foreach (List<SAMAlignedSequence> seqList in readsNotInPloidyForIndividualsDict.Values)
            {
                qualsOut[i++] = seqList[0].MapQ;
            }


            return new double[] 
            { 
                quals.Length > 0 ? Math.Round(quals.Average(), 2) : 0, 
                qualsIn.Length > 0 ? Math.Round(qualsIn.Average(), 2) : 0, 
                qualsOut.Length > 0 ? Math.Round(qualsOut.Average(), 2) : 0
            };
        }

        private double ploidyDisagreement;
        /// <summary>
        /// get the average read quality for all input reads (not distinct reads - every read)
        /// this looks fairly accurate
        /// the qualities are just values, i don't have an upper or lower limit on what they *should* be
        /// </summary>
        /// <param name="dict"></param>
        /// <returns></returns>
        private double[] FindReadQualities()
        {
            double[] readQualScoresAll = new double[CountAll];
            double[] readQualScoresIn = new double[readsInPloidyForIndividuals.Count];
            double[] readQualScoresOut = new double[readsNotInPloidyForIndividuals.Count];


            int i = 0;
            QualitativeSequence qSeq;
            foreach (SAMAlignedSequence seq in sequences)
            {
                qSeq = new QualitativeSequence(SAMDnaAlphabet.Instance, FastQFormatType.Sanger, GetSequence(seq), GetReadQuality(seq));
                readQualScoresAll[i++] = qSeq.GetQualityScores().Average();
                //readQualScoresAll[i++] = GetReadQuality(seq).ToCharArray().Select(x => (int)x - 33).Average(); // [2,4,4,4,5,4,4,4] per sequence
            }

            i = 0;
            foreach (SAMAlignedSequence seq in readsInPloidyForIndividuals)
            {
                qSeq = new QualitativeSequence(SAMDnaAlphabet.Instance, FastQFormatType.Sanger, GetSequence(seq), GetReadQuality(seq));
                readQualScoresIn[i++] = qSeq.GetQualityScores().Average();
                //readQualScoresAll[i++] = GetReadQuality(seq).ToCharArray().Select(x => (int)x - 33).Average(); // [2,4,4,4,5,4,4,4] per sequence
            }

            i = 0;
            foreach (SAMAlignedSequence seq in readsNotInPloidyForIndividuals)
            {
                qSeq = new QualitativeSequence(SAMDnaAlphabet.Instance, FastQFormatType.Sanger, GetSequence(seq), GetReadQuality(seq));
                readQualScoresOut[i++] = qSeq.GetQualityScores().Average();
                //readQualScoresAll[i++] = GetReadQuality(seq).ToCharArray().Select(x => (int)x - 33).Average(); // [2,4,4,4,5,4,4,4] per sequence
            }

           


            return new double[] 
            { 
                readQualScoresAll.Length > 0 ? Math.Round(readQualScoresAll.Average(), 2) : 0, 
                readQualScoresIn.Length > 0 ? Math.Round(readQualScoresIn.Average(), 2) : 0, 
                readQualScoresOut.Length > 0 ? Math.Round(readQualScoresOut.Average(), 2) : 0
            };
        }


        /*private IEnumerable<SAMAlignedSequence> test()
        {
            BAMParser parser = new BAMParser(handler, false);
            parser

            var bamparser = new BAMParser();
            header = bamparser.GetHeader(stream);
            while (!bamparser.IsEOF())
            {
                yield return bamparser.GetAlignedSequence(false);
            }
        }*/

        // Fx with one seq per individual
        private List<int> GetFrequencyDistribution(Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> dict)
        {
            List<int> frequencies = new List<int>();
            foreach (KeyValuePair<String, Dictionary<String, List<SAMAlignedSequence>>> seq in dict)
            {
                frequencies.Add(seq.Value.Count); // where count is the number of sub-dict items
            }
            return frequencies;
        }


        #endregion

        private void IterateSampleDict()
        {
            Console.Write("1");
            int i = 0;
            /*foreach (KeyValuePair<string, List<SAMAlignedSequence>> kvp in sampleDict)
            {


            }*/
            // are all individuals represented equally - for this I need the variation not hte average
            // this is the total read counts per person
            sampleReadCountsAll = new double[CountSamples];
            i = 0;
            foreach (List<SAMAlignedSequence> seqList in sampleDict.Values)
            {
                // for each sample, the total number of reads it has
                sampleReadCountsAll[i++] = seqList.Count;
            }
            Console.Write("2");
            frequencyDistributionSamples = GetFrequencyDistribution(SampleDict);
            Console.Write("3");
        }
        private void IterateSequenceSampleDict()
        {

        }
        private void IterateSampleSequenceDict()
        {
            Console.Write("A");
            clustSeqFrequencies = SampleFrequenciesAvg(sampleSequenceDict);
            Console.Write("B");
            Debug.Assert(Math.Round(clustSeqFrequencies.Sum(), 2) == 1);
            Console.Write("C");

            sampleReadCountsDistinct = new double[CountSamples];
            int i = 0;
            foreach (Dictionary<string, List<SAMAlignedSequence>> seqList in sampleSequenceDict.Values)
            {
                // for each sample, the total number of reads it has
                sampleReadCountsDistinct[i++] = seqList.Count;
            }
            Console.Write("D");
        }
        /// <summary>
        /// Reset all values to null and lists to empty.
        /// </summary>
        public void Reset()
        {
            //sequenceMap.Clear();
            sequences.Clear();
            sequenceDict.Clear();
            sampleDict.Clear();
            sampleSequenceDict.Clear();
            sequenceSampleDict.Clear();
            readsInPloidyForIndividuals.Clear();
            readsNotInPloidyForIndividuals.Clear();

            readsInPloidyForIndividualsDict.Clear();
            readsNotInPloidyForIndividualsDict.Clear();
        }

        #endregion

        #region Private Methods


        /// <summary>
        /// get the number of base pair differences between 2 query strings
        /// loaded with assumptions
        /// </summary>
        /// <param name="seqA"></param>
        /// <param name="seqB"></param>
        /// <returns></returns>
        private int NumQueryDiff(string seqA, string seqB)
        {
            // if diff lengths or diff index, do something about that todo aw
            //int startPos = (seqA.MPos > seqB.MPos) ? seqA.MPos - seqB.MPos : seqB.MPos - seqA.MPos; // startPos is 0 based

            char[] charsA = seqA.ToCharArray();
            char[] charsB = seqB.ToCharArray();
            int count = 0;
            for (int i = 0; i < charsA.Length; i++)
            {
                if (charsA[i] != charsB[i])
                {
                    ++count;
                }
            }
            return count;
        }


        private double GetCountInPloidy(int ploidy, double[] frequencies)
        {
            int count = 0;
            double seqsCount = 0;
            foreach (double fx in frequencies)
            {
                if (count++ < ploidy)
                {
                    seqsCount += fx;
                }
                else
                {
                    break;
                }
            }
            return seqsCount;
        }

        private int GetCountInPloidy(int ploidy, List<int> frequencies)
        {
            int count = 0, seqsCount = 0;
            foreach (int fx in frequencies)
            {
                if (count++ < ploidy)
                {
                    seqsCount += fx;
                }
                else
                {
                    break;
                }
            }
            return seqsCount;
        }

        private double GetPercentBeyondPloidy(int ploidy, List<int> frequencies)
        {
            return Math.Round((frequencies.Sum() - GetCountInPloidy(ploidy, frequencies)) / (Double)frequencies.Sum(), 2);
        }
        

        // is this it??????
        //http://bioinformatics.oxfordjournals.org/content/early/2012/10/09/bioinformatics.bts601.abstract
        //http://bioinformatics.tudelft.nl/

        // gene prediction: uses hmm
        // "copy number variation"

        // "De novo detection of copy number variation "

        // Aneuploidy - abnormal num chromosomes
        // Given that we know the probable ploidy of the organism, we could be said to be
        // detecting aneuoploidy instead

        // each indiv should ONLY have two dominant sequences.
        // for a group of individuals, snps could make it appear that there are >2 dominant sequences


        // {900, 850, 67, 53, 2, 2, 2, 2, 1, 1, 1, 1, 1}
        // avg frequency haploid (top 1)  V avg frequency all the others
        // avg frequency diploid (top 2)
        // avg frequency triploid (top 3)

        // surely we would also expect one or two sequences to be highly represented,
        // and the rest to be lowly represented for true diploid cluster (i.e. bigger gap
        // between 0:1 and 2:max frequency of occurrence)

        // hmm
        // states: haploid, diploid, polyploid, ... (1, 2, 3, 4, more)
        // observations: AAA, AAB, ABA  (i.e. num unique highly represented sequences)
        // start probability: haploid 0.1, diploid 0.7, polyploid 0.2  ??
        // transition: ?????
        // emission probability:    haploid {AAA, AAA 0.9; AAA, AAB 0.1}
        //                          diploid {AAA, AAB 0.8; AAA, AAA 0.1, AAA, AAB, ABA 0.1}

        // trying to observe the actual state based on the observations


        // infer the nmber of haplotypes present within each cluster 
        // (where the reads are the observed data)
        private void HiddenMarkov()
        {
            // don't understand how this translates to transition over time
            // http://nar.oxfordjournals.org/content/33/suppl_2/W451.long need the species name + our seq too short
        }

        // use information derived from the alignment within a cluster as covariates
        private void LogisticRegression()
        {
            // appears to be a sensible model to use
            //(multinomial logitic progression)
        }


        // Get sequence string
        private String GetSequenceDictKey(SAMAlignedSequence seq)
        {
            return GetSequence(seq);
        }

        private string GetSequence(SAMAlignedSequence seq)
        {
            String seqStr = seq.QuerySequence.ToString();
            return Regex.Split(seqStr, "\r\n")[0]; // todo aw something more efficient than regex?
        }

        private string GetReadQuality(SAMAlignedSequence seq)
        {
            String seqStr = seq.QuerySequence.ToString();
            return Regex.Split(seqStr, "\r\n")[1]; // todo aw something more efficient than regex?
        }

        // Get individual string
        private String GetSampleDictKey(SAMAlignedSequence seq)
        {
            foreach (SAMOptionalField field in seq.OptionalFields)
            {
                // I iterate through to find RG each time in case the optional fields
                // do not have a consistent format. Can I assume RG is always OptionalFields[0]
                // to avoid creating iterator and running this loop? todo
                if (field.Tag == "RG")
                {
                    return field.Value;
                }
            }
            return null;
        }

        /// <summary>
        /// Create a sequence dictionary from the current list of aligned sequences
        /// This is a dictionary where each distinct query sequence is represented as a key, and 
        /// a list of all individual sequences that share the same query string is stored
        /// as the value
        /// Each value in the dictionary is a shallow copy of an aligned sequence in sequences
        /// </summary>
        private Dictionary<String, List<SAMAlignedSequence>> MakeSequenceDict(List<SAMAlignedSequence> seqs)
        {
            Dictionary<String, List<SAMAlignedSequence>> dict = new Dictionary<String, List<SAMAlignedSequence>>();

            foreach (SAMAlignedSequence seq in seqs)
            {
                AddToDict(dict, GetSequenceDictKey(seq), seq);
            }
            var sortedSequenceDict = (from sequence in dict
                                      orderby sequence.Value.Count descending
                                      select sequence)
                    .ToDictionary(pair => pair.Key, pair => pair.Value);

            dict = sortedSequenceDict;
            return dict;
        }

        /// <summary>
        /// Create a sample dictionary from the current list of aligned sequences
        /// This is a dictionary where each sample individual is represented as a key, and 
        /// a list of all individual sequences for that individual is stored
        /// as the value
        /// Each value in the dictionary is a shallow copy of an aligned sequence in sequences
        /// </summary>
        private Dictionary<String, List<SAMAlignedSequence>> MakeSampleDict(List<SAMAlignedSequence> seqs)
        {
            Dictionary<String, List<SAMAlignedSequence>> dict = new Dictionary<String, List<SAMAlignedSequence>>();

            foreach (SAMAlignedSequence seq in seqs)
            {
                AddToDict(dict, GetSampleDictKey(seq), seq);
            }
            return dict;
        }

        // Add sequence to dictionary using key
        private void AddToDict(Dictionary<String, List<SAMAlignedSequence>> dict, string key, SAMAlignedSequence seq)
        {
            if (key != null)
            {
                if (dict.ContainsKey(key))
                {
                    List<SAMAlignedSequence> existingVal = dict[key];
                    existingVal.Add(seq);
                    dict[key] = existingVal;
                }
                else
                {
                    dict.Add(key, new List<SAMAlignedSequence> { seq });
                }
            }
            else
            {
                throw new ArgumentException("Invalid dictionary key");
            }
        }

        private Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> MakeNestedSampleDict
            (Dictionary<String, List<SAMAlignedSequence>> fromDict)
        {
            Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> dict
                = new Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>>();

            foreach (KeyValuePair<String, List<SAMAlignedSequence>> entry in fromDict)
            {
                dict[entry.Key] = MakeSampleDict(entry.Value);
            }
            return dict;
        }

        /// <summary>
        /// returns sampleSequenceDict
        /// </summary>
        /// <param name="subDict"></param>
        /// <returns></returns>
        private Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> MakeNestedSequenceDict
            (Dictionary<String, List<SAMAlignedSequence>> subDict)
        {
            readsInPloidyForIndividuals = null;
            readsNotInPloidyForIndividuals = null;

            Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> superDict = new Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>>();

            foreach (KeyValuePair<String, List<SAMAlignedSequence>> entry in subDict) // each sample in sampleDict
            {
                superDict[entry.Key] = MakeSequenceDict(entry.Value); // superDict["sampleName"] = [a sequence dictionary]
                if (readsInPloidyForIndividuals == null)
                {
                    readsInPloidyForIndividuals = new List<SAMAlignedSequence>();
                    readsNotInPloidyForIndividuals = new List<SAMAlignedSequence>();
                }

                int i = 0;
                foreach (List<SAMAlignedSequence> t in superDict[entry.Key].Values) // for each sequence dictionary
                {
                    if(i++ < expectedPloidy){
                        readsInPloidyForIndividuals.AddRange(t);
                    }
                    else
                    {
                        readsNotInPloidyForIndividuals.AddRange(t);
                    }
                }
            }

            readsInPloidyForIndividualsDict = MakeSequenceDict(readsInPloidyForIndividuals);
            readsNotInPloidyForIndividualsDict = MakeSequenceDict(readsNotInPloidyForIndividuals);
            Debug.Assert((readsInPloidyForIndividuals.Count + readsNotInPloidyForIndividuals.Count) == sequences.Count);

            // todo fixme there might be a sequence in both buckets that shares the same query string
            // this is valid but it complicates the averages
            //Console.WriteLine(".." + CountDistinct + "< " + ctDist);

            return superDict;
        }

       


        

        private string PrintSampleFrequencies(Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> sampleDict)
        {
            List<double>[] thisClusterSampleFrequencies = GetSampleFrequenciesAll(sampleDict);
            string returnStr = "";
            int j = 0;
            foreach(List<double> fx in thisClusterSampleFrequencies)
            {
                returnStr += ("Sample "+(j++) + "\t");
                foreach(double d in fx)
                {
                    returnStr += (d + "\t");
                }
                returnStr += ("\n"); // end this list of frequencies for this indiv
                
            }
            return returnStr;
        }

 
        /// <summary>
        /// For an array of values and a given ploidy level, returns an array with three double values, 
        /// representing [averageAllValues, averageInPloidy, averageNotInPloidy]
        /// Values are rounded to 2 decimal places
        /// FIXME is there a copy paste alternative?
        /// </summary>
        /// <param name="arr"></param>
        /// <param name="ploidy"></param>
        /// <returns></returns>
        private double[] PloidyAwareAverages(int[] arr, int ploidy)
        {
            double sumInPloidy = arr.Take(ploidy).Sum();
            double sumOutOfPloidy = arr.Sum() - sumInPloidy;            
            //double avgAll = arr.Sum() / (double)arr.Length;

            double avgAll = (arr.Length > 0) ? arr.Average() : 0;

            double avgIn = sumInPloidy / (double)ploidy;
            double avgOut = (arr.Sum() - sumInPloidy) / (double)(arr.Length - ploidy);


            return new double[] { Math.Round(avgAll, 2), Math.Round(avgIn, 2), Math.Round(avgOut, 2) };


        }
        private double[] PloidyAwareAverages(double[] arr, int ploidy)
        {
            double sumInPloidy = arr.Take(ploidy).Sum();
            double sumOutOfPloidy = arr.Sum() - sumInPloidy;
            double avgAll = (arr.Length > 0) ? arr.Average() : 0;

            double avgIn = sumInPloidy / (double)ploidy;
            double avgOut = (arr.Sum() - sumInPloidy) / (double)(arr.Length - ploidy);
            return new double[] { Math.Round(avgAll, 2), Math.Round(avgIn, 2), Math.Round(avgOut, 2) };
        }


        /// <summary>
        /// Quality score between 0 and 93
        /// </summary>
        /// <param name="readQualStr"></param>
        /// <returns></returns>
        private int[] GetReadQualityScores(string readQualStr)
        {
            int[] scores = new int[readQualStr.Length];
            int i = 0;
            foreach(char c in readQualStr)
            {
                scores[i++] = (int)c - 33;
            }
            return scores;
        }

        

        private double[] SampleFrequenciesAvg(Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> sampleDict)
        {
            
            List<double>[] thisClusterSampleFrequencies = GetSampleFrequenciesAll(sampleDict);
            //Console.WriteLine("The length of the array is " + thisClusterSampleFrequencies.Length);

            int maxlen = 0;
            foreach (List<double> fx in thisClusterSampleFrequencies)
            {
                maxlen = (fx.Count > maxlen) ? fx.Count : maxlen;
            }
            //Console.WriteLine("maxlen: "+maxlen);
            double[] temp = new double[maxlen];


            Array.Clear(temp, 0, maxlen);
            foreach (List<double> fx in thisClusterSampleFrequencies)
            {
                int i;
                for (i = 0; i < fx.Count; i++ )
                {
                    temp[i] += fx[i];
                }
                
            }
            //Console.Write("FX....\t");
            for (int j = 0; j < temp.Length; j++)
            {
                
                temp[j] = temp[j] / (double)thisClusterSampleFrequencies.Length; // divide by num samples
                //Console.Write(temp[j] + "\t");
            }
            //Console.WriteLine("\n");
            return temp;
        }

        private string PoissonifySampleFrequences(List<int>[] frequencies)
        {
            PoissonDistribution[] distributions = new PoissonDistribution[frequencies.Length];
            int i = 0;
            foreach(List<int> sample in frequencies)
            {
                distributions[i++] = new PoissonDistribution(lambda: sample.Average());
            }

            Mixture<PoissonDistribution> mix = new Mixture<PoissonDistribution>(distributions);
            
            double mean = mix.Mean;     // 3.5 
            double median = mix.Median;   // 3.4999998506015895 
            double var = mix.Variance; // 3.25 
            
            double pmf1 = mix.ProbabilityDensityFunction(x: 1);
            double pmf2 = mix.ProbabilityDensityFunction(x: 2);
            double pmf3 = mix.ProbabilityDensityFunction(x: 3);
            double pmf4 = mix.ProbabilityDensityFunction(x: 4);
            double pmf5 = mix.ProbabilityDensityFunction(x: 5);
            double pmf6 = mix.ProbabilityDensityFunction(x: 6);
            return "pmfs: \t" + pmf1 + "\t" + pmf2 + "\t" + pmf3 + "\t" + pmf4 + "\t" + pmf5 + "\t" + pmf6 + "\n";

        }


        private void LogisticRegressionIfy()
        {
            // Suppose we have the following data about some patients. 
            // The first variable is continuous and represent patient 
            // age. The second variable is dicotomic and give whether 
            // they smoke or not (This is completely fictional data). 
            double[][] input =
                {
                    new double[] { 55, 0 }, // 0 - no cancer 
                    new double[] { 28, 0 }, // 0 
                    new double[] { 65, 1 }, // 0 
                    new double[] { 46, 0 }, // 1 - have cancer 
                    new double[] { 86, 1 }, // 1 
                    new double[] { 56, 1 }, // 1 
                    new double[] { 85, 0 }, // 0 
                    new double[] { 33, 0 }, // 0 
                    new double[] { 21, 1 }, // 0 
                    new double[] { 42, 1 }, // 1
                };

            // We also know if they have had lung cancer or not, and  
            // we would like to know whether smoking has any connection 
            // with lung cancer (This is completely fictional data). 
            double[] output =
                {
                    0, 0, 0, 1, 1, 1, 0, 0, 0, 1
                };


            // To verify this hypothesis, we are going to create a logistic 
            // regression model for those two inputs (age and smoking).
            LogisticRegression regression = new LogisticRegression(inputs: 2);

            // Next, we are going to estimate this model. For this, we 
            // will use the Iteravely reweighted least squares method. 
            var teacher = new IterativeReweightedLeastSquares(regression);

            // Now, we will iteratively estimate our model. The Run method returns 
            // the maximum relative change in the model parameters and we will use 
            // it as the convergence criteria. 

            double delta = 0;
            do
            {
                // Perform an iteration
                delta = teacher.Run(input, output);

            } while (delta > 0.001);

            // At this point, we can compute the odds ratio of our variables. 
            // In the model, the variable at 0 is always the intercept term,  
            // with the other following in the sequence. Index 1 is the age 
            // and index 2 is whether the patient smokes or not. 

            // For the age variable, we have that individuals with 
            //   higher age have 1.021 greater odds of getting lung 
            //   cancer controlling for cigarrete smoking. 
            double ageOdds = regression.GetOddsRatio(1); // 1.0208597028836701 

            // For the smoking/non smoking category variable, however, we 
            //   have that individuals who smoke have 5.858 greater odds 
            //   of developing lung cancer compared to those who do not  
            //   smoke, controlling for age (remember, this is completely 
            //   fictional and for demonstration purposes only). 
            double smokeOdds = regression.GetOddsRatio(2); // 5.8584748789881331
        }

        // Gets an array of List<double> types
        // the length of the array is the number of individuals in the cluster
        private List<double>[] GetSampleFrequenciesAll(Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> sampleDict)
        {
            List<double>[] lists = new List<double>[sampleDict.Count];
            int i = 0;
            foreach (Dictionary<String, List<SAMAlignedSequence>> sample in sampleDict.Values) // for each indiv
            {
                lists[i++] = GetSampleFrequencies(sample);
            }
            return lists;
        }

        
        /// <summary>
        /// For a single sample in a single cluster, get top1, top2, ...topAll
        /// </summary>
        /// <param name="dict"></param>
        /// <returns></returns>
        private List<double> GetSampleFrequencies(Dictionary<String, List<SAMAlignedSequence>> dict)
        {
            List<double> frequencies = new List<double>();

            int totalSeqs = 0;
            foreach (List<SAMAlignedSequence> seq in dict.Values)
            {
                totalSeqs += seq.Count;
                
            }

            foreach (List<SAMAlignedSequence> seq in dict.Values)
            {
                frequencies.Add(seq.Count / (double)totalSeqs);
            }

            return frequencies.OrderByDescending(item => item).ToList();
        }




        #endregion


        public void Shrink()
        {
            Console.Write("1");
            sequenceDict.Clear();
            Console.Write("2");
            sampleDict.Clear();
            Console.Write("3");
            sequenceSampleDict.Clear();
            Console.Write("4");
            sampleSequenceDict.Clear();
            Console.Write("5");
            readsInPloidyForIndividuals.Clear();
            Console.Write("6");
            readsNotInPloidyForIndividuals.Clear();
            Console.Write("7");
            
        }

    }
    

    


}
