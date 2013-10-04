using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace Bio.Algorithms.Metric
{
    /// <summary>
    /// An IMetric calculates various metric values from a list of SAMAlignedSequences, where each
    /// list of sequences belongs to a cluster. Clusters are determined based on sequence similarity, 
    /// e.g. as calculated by an MCL graph clustering algorithm, or any other clustering approach
    /// </summary>
    public class ClusterMetric : IMetric
    {

        #region Private Fields

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

        #endregion

        #region Constructors

        /// <summary>
        /// The default constructor.
        /// </summary>
        public ClusterMetric()
        {
            //throw new NotImplementedException();
        }

        #endregion

        #region Properties

        /// <summary>
        /// Cluster ID (unique ID of reference sequence against which all sequences in the cluster
        /// are aligned)
        /// (If sequences currently being handled are from multiple clusters, only the first Id
        /// will be returned)
        /// </summary>
        public string Id
        {
            get 
            {
                if (! sequences[0].Flag.HasFlag(SAMFlags.UnmappedQuery))
                {
                    String rname = sequences[0].RName;
                    return (sequences != null) ? rname : null; 
                }
                else
                {
                    throw new Exception("Unmapped query exception"); // todo particular kind of exception?
                }
                
            }
        }

        /// <summary>
        /// Total number of all sequences in the cluster
        /// </summary>
        public int CountAll
        {
            get { return (sequences != null) ? sequences.Count : 0; }
        }

        /// <summary>
        /// Total number of distinct sequences in the cluster
        /// </summary>
        public int CountDistinct
        {
            get { return (sequenceDict != null) ? sequenceDict.Count : 0; }
        }

        /// <summary>
        /// Total number of samples (individuals) represented in the cluster
        /// </summary>
        public int CountSamples
        {
            get { return (sampleDict != null) ? sampleDict.Count : 0; }
        }

        /// <summary>
        /// A dictionary where each distinct query sequence is represented as a key, and 
        /// a list of all individual sequences that share the same query string is stored
        /// as the value
        /// </summary>
        public Dictionary<String, List<SAMAlignedSequence>> SequenceDict
        {
            get { return sequenceDict; }
        }

        /// <summary>
        /// A dictionary where each sample individual is represented as a key, and 
        /// a list of all individual sequences for that individual is stored
        /// as the value
        /// </summary>
        public Dictionary<String, List<SAMAlignedSequence>> SampleDict
        {
            get { return sampleDict; }
        }

        /// <summary>
        /// List of all sequences.
        /// </summary>
        public List<SAMAlignedSequence> Sequences
        {
            get { return sequences; }
        }
        
        /// <summary>
        /// Experimental.
        /// </summary>
        public List<int> FrequencyDistributionSequences
        {
            get { return GetFrequencyDistribution(SequenceDict); }
        }

        /// <summary>
        /// Experimental.
        /// </summary>
        public List<int> FrequencyDistributionSamples
        {
            get { return GetFrequencyDistribution(SampleDict); }
        }

        /// <summary>
        /// The total number of sequences which are within ploidy (i.e. the top two represented query
        /// sequences)
        /// </summary>
        public int CountInPloidy
        {
            get { return GetCountInPloidy(2, GetFrequencyDistribution(SequenceDict)); }
        }

        /// <summary>
        /// The total number of sequences which are beyond ploidy (i.e. total count excluding the 
        /// top two represented query sequences)
        /// </summary>
        public Double BeyondPloidy
        {
            get { return GetPercentBeyondPloidy(2, GetFrequencyDistribution(SequenceDict)); }
        }

        
        #endregion

        #region Public Methods

        /// <summary>
        /// String of tab-separated values for writing to file by MetricFormatter.
        /// todo fixme we want the MetricFormater to decide how to write each line, but
        /// for now until I know what values are being written I will let the Metric handle this.
        /// </summary>
        public string ToFileString()
        {
            double stdDirt = GetPercentBeyondPloidy(2, GetFrequencyDistribution(SequenceDict));
            
            //double indivDirt = GetPercentBeyondPloidy(2, GetFrequencyDistribution(sequenceSampleDict));
            //double anorDirt = (stdDirt + indivDirt) / 2;


            //List<int> scaledFx = new List<int>();
            //List<int> seqFx = GetFrequencyDistribution(SequenceDict);
            //List<int> indivFx = GetFrequencyDistribution(sequenceSampleDict);
            /*for (int i = 0; i < seqFx.Count; i++)
            {
                scaledFx.Add(seqFx[i] * indivFx[i]);
            }*/
            //double scaledDirt = GetPercentBeyondPloidy(2, scaledFx);

            //Dictionary<double, int> baseFreqDict = BaseFrequencies();
            // for each seq in the cluster

            string baseFxString = "";
            List<Dictionary<char, double>> baseFx = BaseFrequencies();
            foreach(Dictionary<char, double> position in baseFx)
            {
                //baseFxString = baseFxString + "(";
                foreach(KeyValuePair<char, double> pos in position)
                {
                    baseFxString = baseFxString + pos.Key + "\t" + pos.Value + "\t";
                }
                baseFxString = baseFxString + "\n ";
            }

            List<int> frequenciesWSnp = FreqCombineSnp();

            double dirtWSnp = GetPercentBeyondPloidy(2, frequenciesWSnp);

            return Id + "\t" +
                CountDistinct + "\t" +
                CountSamples + "\t" +
                stdDirt + ", "+dirtWSnp
                //+ indivDirt + ", " + anorDirt + ", *" + scaledDirt
                //+ Environment.NewLine
                + "\n" + FrequencyAndGaps() + "\n" //+ baseFxString;
                ;
                //+ CompositeDictMetrics(sampleSequenceDict) 
                //+ CompositeDictMetrics(sequenceSampleDict);
                ;
        }

        /// <summary>
        /// Calculate metric values from the given list of sequences.
        /// </summary>
        /// <param name="sequences">Sequences to add.</param>

        public void Calculate(List<SAMAlignedSequence> sequences)
        {
            this.sequences = sequences;
            sequenceDict = MakeSequenceDict(sequences);
            sampleDict = MakeSampleDict(sequences);
            sampleSequenceDict = MakeNestedSequenceDict(sampleDict);
            sequenceSampleDict = MakeNestedSampleDict(sequenceDict);
        }

        /// <summary>
        /// Reset all values to null and lists to empty.
        /// </summary>
        public void Reset()
        {
            sequences.Clear();
            sequenceDict.Clear();
            sampleDict.Clear();
            sampleSequenceDict.Clear();
            sequenceSampleDict.Clear();
        }

        #endregion

        #region Private Methods

        // Just testing out some values
        private string FrequencyDistributionMetrics()
        {
            List<int> frequencies = FrequencyDistributionSequences;
            int countAll = frequencies.Sum(); 
            int length = frequencies.Count();

            int haploidAverage = frequencies[0] / 1;
            int nonHaploidAverage = (countAll - frequencies[0]) / (length - 1);

            int diploidAverage = (frequencies[0] + frequencies[1]) / 2;
            int nonDiploidAverage = (countAll - frequencies[0] - frequencies[1]) / (length - 2);

            int triploidAverage = (frequencies[0] + frequencies[1] + frequencies[2]) / 3;
            int nonTriploidAverage = (countAll - frequencies[0] - frequencies[1] - frequencies[2]) 
                / (length - 3);

            int tetraploidAverage = (frequencies[0] + frequencies[1] + frequencies[2] + frequencies[3]) / 4;
            int nonTetraploidAverage = (countAll - frequencies[0] - frequencies[1] - frequencies[2] - frequencies[3]) 
                / (length - 4);

            // we are looking at the size of the gap between each ploidy level
            int hapDipGap = haploidAverage - diploidAverage;
            int dipTripGap = diploidAverage - triploidAverage;
            int tripQuadGap = triploidAverage - tetraploidAverage;            
            
            /*return frequenciesList + "\n"
                + haploidAverage + "\t" + diploidAverage + "\t" + tetraploidAverage + "\t" + quadraploidAverage + "\n"
                + nonHaploidAverage + "\t" + nonDiploidAverage + "\t" + nonTetraploidAverage + "\t" + nonQuadraploidAverage + "\n"
                + (haploidAverage - nonHaploidAverage) + "\t" + (diploidAverage - nonDiploidAverage) + "\t"
                + (tetraploidAverage - nonTetraploidAverage) + "\t" + (quadraploidAverage - nonQuadraploidAverage) + "\n"
                + "("+hapDipGap + ")\t(" + dipTripGap + ")\t(" + tripQuadGap + ")\n"
                + looksOK + "\n\n"
                ;*/
            return "";
        }


        // loaded with assumptions
        private int NumQueryDiff(string seqA, string seqB)
        {
            // if diff lengths or diff index, do something about that todo aw
            //int startPos = (seqA.MPos > seqB.MPos) ? seqA.MPos - seqB.MPos : seqB.MPos - seqA.MPos; // startPos is 0 based
            
            char[] charsA = seqA.ToCharArray();
            char[] charsB = seqB.ToCharArray();
            int count = 0;
            for (int i = 0; i < charsA.Length; i++ )
            {
                if(charsA[i] != charsB[i])
                {
                    ++count;
                }
            }
            return count;
        }

        private List<int> FreqCombineSnp()
        {
            List<int> frequencies = new List<int>();

            int firstI = 0, secondI = 0;
            int firstSeqIndex = 0, nextSeqIndex;


            for (int i = firstSeqIndex+1; i < SequenceDict.Count; i++)
            {
                int diff = NumQueryDiff(SequenceDict.ElementAt(firstSeqIndex).Key, SequenceDict.ElementAt(i).Key);
                // if 1 snp suggested
                if(diff < 2) 
                {
                    frequencies.Add(SequenceDict.ElementAt(firstSeqIndex).Value.Count + SequenceDict.ElementAt(i).Value.Count);
                    firstI = i;
                    break;
                }
                // worthless sequence or finished iterating
                else if (SequenceDict.ElementAt(i).Value.Count < 2 || i == SequenceDict.Count - 1)
                {
                    frequencies.Add(SequenceDict.ElementAt(firstSeqIndex).Value.Count); // first elem only
                    firstI = -1;
                    break;
                }
            }
            

            nextSeqIndex = (firstI == 1) ? 2 : 1;

            for (int i = nextSeqIndex + 1; i < SequenceDict.Count; i++)
            {
                int diff = NumQueryDiff(SequenceDict.ElementAt(nextSeqIndex).Key, SequenceDict.ElementAt(i).Key);
                // if 1 snp suggested
                if (diff < 2)
                {
                    frequencies.Add(SequenceDict.ElementAt(nextSeqIndex).Value.Count + SequenceDict.ElementAt(i).Value.Count);
                    secondI = i;
                    break;
                }
                // worthless sequence or finished iterating
                else if (SequenceDict.ElementAt(i).Value.Count < 2 || i == SequenceDict.Count - 1)
                {
                    frequencies.Add(SequenceDict.ElementAt(nextSeqIndex).Value.Count); // first elem only
                    secondI = -1;
                    break;
                }
            }

            for (int i = nextSeqIndex + 1; i < SequenceDict.Count; i++)
            {
                if (i != firstI && i != secondI)
                {
                    frequencies.Add(SequenceDict.ElementAt(i).Value.Count);
                }
            }
            return frequencies;
        }


        private string FrequencyAndGaps()
        {
            string returnVal = "";
            
            // indiv fx
            /*List<int> frequencies2 = GetFrequencyDistribution(sequenceSampleDict);
            string frequenciesList2 = string.Join(",", frequencies2.ToArray());           
            returnVal = returnVal + "Indiv fx: "+frequenciesList2 + "\n\n";*/
            
            List<int> frequencies = FrequencyDistributionSequences;
            string frequenciesList = string.Join(",", frequencies.ToArray());

            List<int> frequenciesWSnp = FreqCombineSnp();
            string frequenciesWSnpList = string.Join(",", frequenciesWSnp.ToArray());



            // scaled fx
            /*List<int> frequencies3 = new List<int>();
            for (int i = 0; i < frequencies.Count; i++)
            {
                frequencies3.Add(frequencies2[i] * frequencies[i]);
            }
            string frequenciesList3 = string.Join(",", frequencies3.ToArray());
            returnVal = returnVal + "Scaled fx: " + frequenciesList3 + "\n";*/

            returnVal = returnVal + "Seq fx:\n" + frequenciesList + "\nSeq fx w/ snp:\n" + frequenciesWSnpList + "\n";

            return returnVal;
            /*List<int> gaps = new List<int>();
            
            for(int i = 1; i < frequencies.Count; i++)
            {
                gaps.Add(frequencies[i-1] - frequencies[i]);
                
            }
            string gapsList = string.Join(",", gaps.ToArray());
            int gapCount = gaps.Sum();
            for (int i = 0; i < gaps.Count; i++)
            {
                gaps[i] = gaps[i]*10000/gapCount;
            }
            string normGapsList = string.Join(",", gaps.ToArray());
            return returnVal + " Gaps: " + gapsList + "\n" + normGapsList+"\n";*/
        }


        private string CompositeDictMetricsBySequence()
        {
            return CompositeDictMetrics(sequenceSampleDict);
        }

        private string CompositeDictMetricsBySample()
        {
            return CompositeDictMetrics(sampleSequenceDict);
        }

        private string CompositeDictMetrics(Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> baseDict)
        {
            StringBuilder sb = new StringBuilder();
            // For each item in base dict
            int count = 0;
            foreach (KeyValuePair<String, Dictionary<String, List<SAMAlignedSequence>>> sample in baseDict)
            {
                sb.Append(count.ToString() + ": ");
                foreach (KeyValuePair<String, List<SAMAlignedSequence>> indiv in sample.Value)
                {
                    // the number of SAMAlignmentSequences for each dict-sub-item
                    sb.Append(indiv.Value.Count.ToString()).Append(',');
                }
                sb.AppendLine("");
                count++;
            }
            return sb.ToString();
        }


        // get count in ploidy for distinct sequences
        /*private int GetCountInPloidy(int ploidy)
        {
            int count = 0, seqsCount = 0;
            foreach (KeyValuePair<String, List<SAMAlignedSequence>> seqs in sequenceDict)
            {
                if(count++ < ploidy)
                {
                    seqsCount += seqs.Value.Count;
                }
                else
                {
                    break;
                }
            }
            return seqsCount;
        }*/

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
            Dictionary<String, List<SAMAlignedSequence>>  dict = new Dictionary<String, List<SAMAlignedSequence>>();

            foreach (SAMAlignedSequence seq in seqs)
            {
                AddToDict(dict, GetSequenceDictKey(seq), seq);
            }

            // todo aw inefficient?
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

        private Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> MakeNestedSequenceDict
            (Dictionary<String, List<SAMAlignedSequence>> fromDict)
        {
            Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> dict
                = new Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>>();

            foreach (KeyValuePair<String, List<SAMAlignedSequence>> entry in fromDict)
            {
                dict[entry.Key] = MakeSequenceDict(entry.Value);
            }
            return dict;
        }

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

        // Fx with one seq per individual
        private List<int> GetFrequencyDistribution(Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> dict)
        {
            // todo aw probably a more efficient way
            List<int> frequencies = new List<int>();
            foreach (KeyValuePair<String, Dictionary<String, List<SAMAlignedSequence>>> seq in dict)
            {
                frequencies.Add(seq.Value.Count); // the number of sub-dict items
            }
            return frequencies;
        }



        private List<Dictionary<char, double>> BaseFrequencies()
        {
            //Dictionary<double, int> freqDict = new Dictionary<double, int>();
            List<Dictionary<char, double>> freqList = new List<Dictionary<char, double>>();

            //double A = 0, T = 0, C = 0, G = 0;
            foreach(SAMAlignedSequence seq in Sequences)
            {
                string seqStr = Regex.Split(seq.QuerySequence.ToString(), "\r\n")[0];
                int count = 0;
                foreach (char c in seqStr.ToUpper().ToCharArray())
                {
                    // 5 // ok to index 4
                    if(count >= freqList.Count)
                    {
                        freqList.Add(new Dictionary<char, double>());
                    }
                    if(freqList[count].ContainsKey(c))
                    {
                        ++freqList[count][c];
                    }
                    else
                    {
                        freqList[count].Add(c, 1);
                    }

                    count++;
                }
            }
                
            foreach(Dictionary<char, double> bases in freqList)
            {
                double numBases = bases.Values.Sum();
                if(bases.ContainsKey('A'))
                {
                    bases['A'] = Math.Round((bases['A'] / numBases), 2);
                    if(bases['A'] == 0) { bases.Remove('A'); }
                }
                if (bases.ContainsKey('T'))
                {
                    bases['T'] = Math.Round((bases['T'] / numBases), 2);
                    if (bases['T'] == 0) { bases.Remove('T'); }
                }
                if (bases.ContainsKey('C'))
                {
                    bases['C'] = Math.Round((bases['C'] / numBases), 2);
                    if (bases['C'] == 0) { bases.Remove('C'); }
                }
                if (bases.ContainsKey('G'))
                {
                    bases['G'] = Math.Round((bases['G'] / numBases), 2);
                    if (bases['G'] == 0) { bases.Remove('G'); }
                }
                

                /*var sortedSequenceDict = (from b in bases
                                          orderby b.Key ascending
                                          select b)
                    .ToDictionary(pair => pair.Key, pair => pair.Value);
                bases = sortedSequenceDict;*/

            }

            return freqList;
        }


        #endregion

    }


    
}
