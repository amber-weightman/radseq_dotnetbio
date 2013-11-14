using Accord.Statistics.Distributions.Univariate;
using Accord.Statistics.Models.Regression;
using Accord.Statistics.Models.Regression.Fitting;
using Bio;
using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Text.RegularExpressions;

namespace Ploidulator
{
    /// <summary>
    /// An IMetric calculates various metric values from a list of SAMAlignedSequences, where each
    /// list of sequences belongs to a cluster. Clusters are determined based on sequence similarity, 
    /// e.g. as calculated by an MCL graph clustering algorithm, or any other clustering approach
    /// </summary>
    public class ClusterMetric : IMetric
    {
        #region Private Static Fields
        private static CultureInfo CI = new CultureInfo("en-AU");
        #endregion

        #region Private Fields

        #region general cluster statistics

        /// <summary>
        /// The expected ploidy level of the organism under study (default 2)
        /// </summary>
        private int expectedPloidy = 2;

        /// <summary>
        /// The number of individuals with samples present in the input file
        /// </summary>
        private int numSamples;

        /// <summary>
        /// Cluster ID (unique ID of reference sequence against which all sequences in the cluster
        /// are aligned)
        /// </summary>
        private string id;

        /// <summary>
        /// Cluster reference sequence (against which all sequences in the cluster are aligned)
        /// </summary>
        private string referenceSequence;

        #endregion

        #region sequence storage

        /// <summary>
        /// A dictionary of all reads which are in-ploidy for the individual to which they belong.
        /// Each distinct read sequence is the dictionary key and the value list contains a all copies of that distinct sequence
        /// The value list is sorted in descending order of list size
        /// </summary>
        private Dictionary<string, List<SAMAlignedSequence>> readsInPloidyForIndividualsDict = null;

        /// <summary>
        /// A dictionary of all reads which are NOT in-ploidy for the individual to which they belong.
        /// Each distinct read sequence is the dictionary key and the value list contains a all copies of that distinct sequence
        /// The value list is sorted in descending order of list size
        /// </summary>
        private Dictionary<string, List<SAMAlignedSequence>> readsNotInPloidyForIndividualsDict = null;

        /// <summary>
        /// A dictionary where each distinct query sequence is represented as a key, and 
        /// a list of all  sequences that share the same query string is stored
        /// as the value
        /// The value list is sorted in descending order of list size
        /// </summary>
        private Dictionary<String, List<SAMAlignedSequence>> sequenceDict = null;

        /// <summary>
        /// A dictionary where each sample individual is represented as a key, and 
        /// a list of all  sequences for that individual is stored
        /// as the value.
        /// The value list is sorted in descending order of list size.
        /// </summary>
        private Dictionary<String, List<SAMAlignedSequence>> sampleDict = null;

        /// <summary>
        /// A dictionary which is ordered by individual in the first instance, and is then ordered by distinct query
        /// sequence string, with a copy of each read for that individual/that sequence following.
        /// The value list is sorted in descending order of list size.
        /// </summary>
        private Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> sampleSequenceDict = null;

        /// <summary>
        /// Sequences are initially stored in this collection before they are processed. At the time of processing, will
        /// contain a copy of every read in the cluster
        /// </summary>
        private Collection<SAMAlignedSequence> sequencesThisCluster = null;

        #endregion

        #region phase

        /// <summary>
        /// Formatted string as required as input by PHASE. For each allele position, 'S' indicates biallelic and
        /// 'M' indicates multiallelic 
        /// </summary>
        private string loci;

        /// <summary>
        /// Formatted string as required as input by PHASE. Multi-line string where each individual has three lines: 
        /// unique identifier, first allele at each position, second allele at each position (genotype)
        /// </summary>
        private string phaseData;

        #endregion

        #region stats lists

        /// <summary>
        /// Ordered list of the frequencies for each sequence, from most to least frequent (sums to 1). 
        /// Ideally top [expectedPloidy] sequences would together be ~1 (100%)
        /// </summary>
        private Collection<double> individualSequenceDistributions = null;

        #endregion

        #region counts

        /// <summary>
        /// Number of sequences in the cluster.
        /// </summary>
        private int countAll;

        /// <summary>
        /// Number of distinct sequences in the cluster
        /// </summary>
        private int countDistinct;

        /// <summary>
        /// Number of samples (individuals) represented in the cluster.
        /// </summary>
        private int countSamples;

        /// <summary>
        /// The total read count for each individual.
        /// </summary>
        private Collection<double> sampleReadCountsAll;

        /// <summary>
        /// The number of distinct read counts for each individual.
        /// </summary>
        private Collection<double> sampleReadCountsDistinct;

        #endregion

        #region stats values

        /// <summary>
        /// The number of haplotypes in this cluster (value of -1 indicates haplotypes were not calculated)
        /// </summary>
        private int numberOfHaplotypes = -1;

        /// <summary>
        /// Number of distinct sequences which are in ploidy for some individuals but out of ploidy for others (for the same sequence).
        /// </summary>
        private double ploidyDisagreement = -1;

        /// <summary>
        /// Identifies whether the cluster is good or bad (note that this value can be freely get and set - external class
        /// ClusterMetricHandler is responsible for determining what makes a cluster 'good')
        /// </summary>
        private bool good = false;

        /// <summary>
        /// Cluster 'dirt' (proportion of reads which are outside [expected ploidy] for their individual
        /// </summary>
        private double dirt = -1;

        /// <summary>
        /// Percentage of the entire population which is represented in this cluster. (Precondition: numSamples must
        /// have been set correctly via the constructor)
        /// </summary>
        private double populationPercentage = -1;

        /// <summary>
        /// The average alignment quaity for distinct reads in the cluster.
        /// </summary>
        private double alignmentQuality = -1;

        /// <summary>
        /// The read alignment quaity for distinct reads in the cluster.
        /// </summary>
        private double readQuality = -1;

        #endregion

        #endregion

        #region Properties

        /// <summary>
        /// Read only. Gets the Cluster ID (unique ID of reference sequence against which all sequences in the cluster
        /// are aligned).
        /// </summary>
        public string Id { get { return id; } }

        #region phase formatted strings

        /// <summary>
        /// Read only. Gets the formatted string as required as input by PHASE. For each allele position, 'S' indicates biallelic and
        /// 'M' indicates multiallelic.
        /// </summary>
        public string PhaseLoci { get { return loci; } }

        /// <summary>
        /// Read only. Gets the formatted string as required as input by PHASE. Multi-line string where each individual has three lines: 
        /// unique identifier, first allele at each position, second allele at each position (genotype).
        /// </summary>
        public string PhaseData { get { return phaseData; } }

        #endregion

        #region counts

        /// <summary>
        /// Read only. Gets the number of sequences in the cluster.
        /// </summary>
        public int CountAll { get { return countAll; } }

        /// <summary>
        /// Read only. Gets the number of distinct sequences in the cluster.
        /// </summary>
        public int CountDistinct { get { return countDistinct; } }

        /// <summary>
        /// Read only. Gets the number of samples (individuals) represented in the cluster.
        /// </summary>
        public int CountSamples { get { return countSamples; } }

        #endregion

        #region sequence storage

        /// <summary>
        /// Gets a dictionary where each distinct query sequence is represented as a key, and 
        /// a list of all individual sequences that share the same query string is stored
        /// as the value. Read only.
        /// </summary>
        public Dictionary<String, List<SAMAlignedSequence>> SequenceDictionary { get { return sequenceDict; } }

        /// <summary>
        /// Gets a dictionary where each sample individual is represented as a key, and 
        /// a list of all  sequences for that individual is stored
        /// as the value.
        /// The value list is sorted in descending order of list size. Read only.
        /// </summary>
        public Dictionary<String, List<SAMAlignedSequence>> SampleDictionary { get { return sampleDict; } }

        /// <summary>
        /// Gets a list of all sequences. Read only.
        /// </summary>
        public Collection<SAMAlignedSequence> Sequences { get { return sequencesThisCluster; } }

        #endregion

        #region stats lists

        /// <summary>
        /// Gets an ordered list of the frequencies for each sequence, from most to least frequent (sums to 1). 
        /// Ideally top [expectedPloidy] sequences would together be ~1 (100%). Read only.
        /// </summary>
        public Collection<double> ClusterSequenceFrequencies { get { return individualSequenceDistributions; } }

        /// <summary>
        /// Gets the total read count for each individual. Read only.
        /// </summary>
        public Collection<double> SampleReadCountsAll { get { return sampleReadCountsAll; } }

        /// <summary>
        /// Gets the number of distinct read counts for each individual. Read only.
        /// </summary>
        public Collection<double> SampleReadCountsDistinct { get { return sampleReadCountsDistinct; } }

        #endregion

        #region stats values

        /// <summary>
        /// Gets the number of probable haplotypes in this cluster (value of -1 indicates haplotypes were not calculated). Read only.
        /// </summary>
        public int NumberOfHaplotypes { get { return numberOfHaplotypes; } set { numberOfHaplotypes = value; } }

        /// <summary>
        /// Identifies whether the cluster is good or bad (note that this value can be freely get and set - external class
        /// ClusterMetricHandler is responsible for determining what makes a cluster 'good'). False by default. Read only.
        /// </summary>
        public bool Good { get { return good; } set { good = value; } }

        /// <summary>
        /// Gets the cluster 'dirt' (proportion of reads which are outside [expected ploidy] for their individual. Read only.
        /// </summary>
        public double Dirt { get { return dirt; } }

        /// <summary>
        /// Gets the average alignment quaity for distinct reads in the cluster. Read only.
        /// </summary>
        public double AlignmentQuality { get { return alignmentQuality; } }

        /// <summary>
        /// Gets the percentage of the entire population which is represented in this cluster. (Precondition: numSamples must
        /// have been set to the correct value via the constructor). Read only.
        /// </summary>
        public double PopulationPercentage { get { return populationPercentage; } }

        /// <summary>
        /// Gets the read alignment quaity for distinct reads in the cluster. Read only.
        /// </summary>
        public double ReadQuality { get { return readQuality; } }

        #endregion

        #endregion

        #region Constructors

        /// <summary>
        /// The default constructor.
        /// </summary>
        public ClusterMetric()
        {
        }

        /// <summary>
        /// Non-default constructor. Please use this one.
        /// </summary>
        /// <param name="ploidy">Expected ploidy level of the organism.</param>
        /// <param name="samples">Number of samples in data file.</param>
        public ClusterMetric(int ploidy, int samples)
        {
            expectedPloidy = ploidy;
            numSamples = samples;
        }

        #endregion

        #region Private Methods

        #region calculate specific metrics

        /// <summary>
        /// Set the value of fields id and referenceSequence
        /// </summary>
        private void SetClusterReferenceIdAndSequence()
        {
            id = GetId(sequencesThisCluster[0]);
            referenceSequence = GetSequence(sequencesThisCluster[0]);
        }

        /// <summary>
        /// Find the percentage of distinct reads which are both in ploidy (top [n]) for some individuals and 
        /// outside ploidy for other individuals (ideally the value should be 0)
        /// </summary>
        private void SetPloidyDisagreement()
        {
            ploidyDisagreement = Math.Round(((readsInPloidyForIndividualsDict.Count
                + readsNotInPloidyForIndividualsDict.Count) - CountDistinct) / (double)CountDistinct, 2);
        }


        /// <summary>
        /// Initialises individualSequenceDistributions, which represents the distribution of each distinct sequence,
        /// for each individual (this distribution would be used to calculate dirt per individual)
        /// </summary>
        private void SetIndividualSequenceDistributions()
        {
            List<double>[] thisClusterSampleFrequencies = new List<double>[sampleSequenceDict.Count];
            int i = 0;
            int maxlen = 0;
            foreach (Dictionary<String, List<SAMAlignedSequence>> sample in sampleSequenceDict.Values) // for each individual
            {
                thisClusterSampleFrequencies[i] = GetSampleFrequencies(sample);
                maxlen = (thisClusterSampleFrequencies[i].Count > maxlen) ? thisClusterSampleFrequencies[i].Count : maxlen;
                i++;
            }

            double[] frequencies = new double[maxlen];
            Array.Clear(frequencies, 0, maxlen); // fill with 0s
            foreach (List<double> fx in thisClusterSampleFrequencies)
            {
                for (int j = 0; j < fx.Count; j++)
                {
                    frequencies[j] += fx[j];
                }
            }
            for (int j = 0; j < frequencies.Length; j++)
            {
                frequencies[j] = frequencies[j] / (double)thisClusterSampleFrequencies.Length;
            }

            // Convert frequencies to a collection (type more commonly used in this class)
            Collection<double> frequenciesAsCollection = new Collection<double>();
            foreach (double f in frequencies)
            {
                frequenciesAsCollection.Add(f);
            }
            individualSequenceDistributions = frequenciesAsCollection;
        }

        #endregion


        #region iterate dictionaries

        /// <summary>
        /// Iterate through sequenceDictionary once only, and perform operations on each sequence list
        /// Set the average read and alignment qualities for all distinct reads, and the frequency of occurrence of each
        /// distinct read
        /// </summary>
        private void IterateSequenceDict()
        {
            double[] alignmentQualities = new double[CountDistinct]; // one for every sequence in the map, just get its qualities
            double[] readQualities = new double[CountDistinct]; //
            Collection<int> frequencies = new Collection<int>();

            int i = 0;
            foreach (List<SAMAlignedSequence> seqList in sequenceDict.Values)
            {
                // Alighment qualities
                alignmentQualities[i] = seqList[0].MapQ;

                // Read qualities
                QualitativeSequence qSeq = new QualitativeSequence(SAMDnaAlphabet.Instance, FastQFormatType.Sanger, GetSequence(seqList[0]), GetReadQuality(seqList[0]));
                readQualities[i++] = qSeq.GetQualityScores().Average();

                // Frequencies
                frequencies.Add(seqList.Count);
            }

            alignmentQuality = alignmentQualities.Length > 0 ? Math.Round(alignmentQualities.Average(), 2) : 0;
            readQuality = readQualities.Length > 0 ? Math.Round(readQualities.Average(), 2) : 0;

            //frequencyDistributionSequences = frequencies;
        }



        /// <summary>
        /// Iterate through sampleDictionary once only, and perform operations on each sequence list
        /// </summary>
        private void IterateSampleDict()
        {
            // For each individual, for each distinct sequence they have, get the count of exactly how many
            // of each distinct sequence they have, and store these counts in sampleReadCountsAll
            sampleReadCountsAll = new Collection<double>();
            foreach (List<SAMAlignedSequence> seqList in sampleDict.Values)
            {
                sampleReadCountsAll.Add(seqList.Count);
            }
        }


        /// <summary>
        /// Sets cluster dirt. Also produces a count of the number of distinct reads each individual has
        /// </summary>
        private void IterateSampleSequenceDict()
        {
            SetIndividualSequenceDistributions();
            dirt = Math.Round(1 - GetCountInPloidy(expectedPloidy, individualSequenceDistributions), 2);
            // dirt = readsNotInPloidyForIndividualsDict.Count / (double)(readsInPloidyForIndividualsDict.Count + readsNotInPloidyForIndividualsDict.Count);

            sampleReadCountsDistinct = new Collection<double>();
            foreach (Dictionary<string, List<SAMAlignedSequence>> seqList in sampleSequenceDict.Values)
            {
                sampleReadCountsDistinct.Add(seqList.Count);
            }
        }

        #endregion

        #region make/initialise dictionaries


        /// <summary>
        /// returns sampleSequenceDict
        /// adds a nested sequence dict to sampleDict
        /// </summary>
        /// <param name="sampleDict"></param>
        /// <returns></returns>
        private void MakeSampleSequenceDict()
        {
            Collection<SAMAlignedSequence> readsInPloidyForIndividuals = null;
            Collection<SAMAlignedSequence> readsNotInPloidyForIndividuals = null;

            sampleSequenceDict = new Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>>();

            foreach (KeyValuePair<String, List<SAMAlignedSequence>> entry in sampleDict) // each sample in sampleDict
            {
                sampleSequenceDict[entry.Key] = MakeSequenceDict(new Collection<SAMAlignedSequence>(entry.Value)); // superDict["sampleName"] = [a sequence dictionary]
                if (readsInPloidyForIndividuals == null)
                {
                    readsInPloidyForIndividuals = new Collection<SAMAlignedSequence>();
                    readsNotInPloidyForIndividuals = new Collection<SAMAlignedSequence>();
                }

                // for each distinct sequence, by individual, add that sequence to one of two categories:
                // either it is in ploidy (top n) for that individual or it is outside ploidy for
                // that individual
                int i = 0;
                foreach (List<SAMAlignedSequence> t in sampleSequenceDict[entry.Key].Values)
                {
                    if (i++ < expectedPloidy)
                    {
                        foreach (SAMAlignedSequence s in t)
                        {
                            readsInPloidyForIndividuals.Add(s);
                        }
                    }
                    else
                    {
                        foreach (SAMAlignedSequence s in t)
                        {
                            readsNotInPloidyForIndividuals.Add(s);
                        }
                    }
                }
            }
            readsInPloidyForIndividualsDict = MakeSequenceDict(readsInPloidyForIndividuals);
            readsNotInPloidyForIndividualsDict = MakeSequenceDict(readsNotInPloidyForIndividuals);
            Debug.Assert((readsInPloidyForIndividuals.Count + readsNotInPloidyForIndividuals.Count) == sequencesThisCluster.Count);
            SetPloidyDisagreement();
        }


        /// <summary>
        /// Create a sample dictionary from the current list of aligned sequences
        /// This is a dictionary where each sample individual is represented as a key, and 
        /// a list of all individual sequences for that individual is stored
        /// as the value
        /// </summary>
        private void MakeSampleDict()
        {
            sampleDict = new Dictionary<String, List<SAMAlignedSequence>>();

            foreach (SAMAlignedSequence seq in sequencesThisCluster)
            {
                AddToDict(sampleDict, GetRgTag(seq), seq);
            }
        }

        #endregion



        #region base position allele haplotype stuff todo fixme

        /// <summary>
        /// Given alleleFxThisIndiv representing each allele and the frequency with which it occurrs, returns
        /// just a list of the allele chars (with unrecognised chars coded as '?')
        /// </summary>
        /// <param name="alleleFxThisIndiv"></param>
        /// <param name="locusIndex"></param>
        /// <returns></returns>
        private char[] GetAllelesAtLocusForIndiv(Dictionary<char, double>[] alleleFxThisIndiv, int locusIndex)
        {
            char[] allelesThisIndiv = alleleFxThisIndiv[locusIndex].Keys.ToArray();
            for (int j = 0; j < allelesThisIndiv.Length; j++)
            {
                if (!alleles.ContainsKey(allelesThisIndiv[j]))
                {
                    allelesThisIndiv[j] = '?';
                }
            }
            return allelesThisIndiv;
        }



        // Dictionary of all possible allele values
        private Dictionary<char, string> alleles = new Dictionary<char, string>()
            {
                {'A', "1"},
                {'T', "2"},
                {'C', "3"},
                {'G', "4"},
                {'?', "-1"}
            };

        
        private static bool GetIndivBiAllelicLocusAlleles(char[] allelesThisIndividual, int numSeqsThisIndivHas, ref string chr1, ref string chr2)
        {
            switch (allelesThisIndividual.Length)
            {
                case 0:
                    // This individual has no alleles at this position
                    chr1 += ("? ");
                    chr2 += ("? ");
                    break;

                case 1:
                    if (numSeqsThisIndivHas > 1)
                    {
                        // This individual has one distinct allele at this position, but more than one instance of it occurring,
                        // so we assume it may be on different chromosomes
                        chr1 += allelesThisIndividual[0] + " ";
                        chr2 += allelesThisIndividual[0] + " ";
                    }
                    else
                    {
                        // Indiv only has one sequence in this cluster. Not enough information to infer that chr2 == chr1 at this position
                        chr1 += allelesThisIndividual[0] + " ";
                        chr2 += "? ";
                    }
                    break;

                case 2:
                    // Indiv has two sequences, both with different alleles at this position
                    chr1 += allelesThisIndividual[0] + " ";
                    chr2 += allelesThisIndividual[1] + " ";
                    break;
                
                default:
                    Console.WriteLine(Properties.Resources.ALLELE_COUNT_ERROR);
                    return false;
            }
            return true;
        }


        private bool GetIndivMultiAllelicLocusAlleles(char[] allelesThisIndividual, int numSeqsThisIndivHas, ref string chr1, ref string chr2)
        {
            switch (allelesThisIndividual.Length)
            {
                case 0:
                    // This individual has no alleles at this position
                    chr1 += ("-1 ");
                    chr2 += ("-1 ");
                    break;

                case 1:
                    if (numSeqsThisIndivHas > 1)
                    {
                        // This individual has one distinct allele at this position, but more than one instance of it occurring,
                        // so we assume it may be on different chromosomes
                        chr1 += alleles[allelesThisIndividual[0]] + " ";
                        chr2 += alleles[allelesThisIndividual[0]] + " ";
                    }
                    else
                    {
                        // Indiv only has one sequence in this cluster. Not enough information to infer that chr2 == chr1 at this position
                        chr1 += alleles[allelesThisIndividual[0]] + " ";
                        chr2 += "-1 ";
                    }
                    break;

                case 2:
                    // Indiv has two sequences, both with different alleles at this position
                    chr1 += alleles[allelesThisIndividual[0]] + " ";
                    chr2 += alleles[allelesThisIndividual[1]] + " ";
                    break;

                default:
                    Console.WriteLine(Properties.Resources.ALLELE_COUNT_ERROR);
                    return false;
            }
            return true;
        }


        


        

        private void ConstructPhaseSnpString()
        {
            
            Dictionary<char, double>[] alleleFxAllIndiv = new Dictionary<char, double>[referenceSequence.Length];
            List<Dictionary<char, double>[]> alleleFxAllIndivFull = new List<Dictionary<char, double>[]>();

            // Generates the data for each individual for the phase file
            string phaseFileData = "";
            foreach (Dictionary<string, List<SAMAlignedSequence>> seqList in sampleSequenceDict.Values) 
            {
                Dictionary<char, double>[] alleleFxThisIndiv = BaseFrequencies(seqList.Keys.ToArray(), ref alleleFxAllIndiv);
                alleleFxAllIndivFull.Add(alleleFxThisIndiv);
            }

            loci = "";
            string noSnp = "-";
            string biAllelic = "S";
            string multiAllelic = "M";
            foreach (Dictionary<char, double> i in alleleFxAllIndiv)
            {
                if(i.Keys.Count == 0 || i.Keys.Count == 1)
                {
                    loci += noSnp;
                } 
                else if(i.Keys.Count == 2)
                {
                    loci += biAllelic;
                }
                else
                {
                    loci += multiAllelic;
                }
            }

            int j = 0;
            foreach (Dictionary<string, List<SAMAlignedSequence>> seqList in sampleSequenceDict.Values)
            {
                // need to construct another dict for each snp pos which exists, to record s or m

                string indivId = "#" + seqList.Values.ToArray()[0][0].QName;
                string chr1 = "", chr2 = "";
                int locusCount = 0;
                foreach (char allele in loci) // for each allele at this locus, where loci is in the format "S--SS-SM---MSSMMS-MMM--MSSSMS"
                    // we will iterate through looking at one physical locus at a time and examining all sequences that align to that locus
                {
                    if (allele != '-')
                    {

                        // Get an array of all alleles which appear at this locus for this individual
                        // If any of these allele characters are not in the dictionary of recognised alleles, replace them with '?' 
                        char[] allelesThisIndiv = GetAllelesAtLocusForIndiv(alleleFxAllIndivFull[j], locusCount);

                        // Between all samples, there are two alleles at this position
                        if (allele == 'S')
                        {
                            GetIndivBiAllelicLocusAlleles(allelesThisIndiv, seqList.Count, ref chr1, ref chr2);
                        }

                        // Between all samples, there are three or four alleles at this position (although this particular individual
                        // will still only have 2)
                        else if (allele == 'M')
                        {
                            GetIndivMultiAllelicLocusAlleles(allelesThisIndiv, seqList.Count, ref chr1, ref chr2);
                        }
                    }
                    locusCount++;
                }

                // Add this individual's genotype information to the phase file data string
                phaseFileData += (indivId + "\r\n" + chr1 + "\r\n" + chr2 + "\r\n");
                j++;
            }

            phaseData = phaseFileData;

            string locii = "";
            for (int k = 0; k < loci.Length; k++ )
            {
                if(loci[k] != '-')
                {
                    locii += loci[k];
                }
            }
            loci = locii;


            //readQuality = FindReadQuality(); // takes a long time

        }


  

        

        // precondition seqs is an ordered list
        /// <summary>
        /// 
        /// </summary>
        /// <param name="seqs"></param>
        /// <param name="masterFreqList"></param>
        /// <returns></returns>
        private Dictionary<char, double>[] BaseFrequencies(string[] seqs, ref Dictionary<char, double>[] masterFreqList)
        {
            
            // A dictionary of char:double for each base pair position
            Dictionary<char, double>[] freqList = new Dictionary<char, double>[seqs[0].Length];

            //For each position, get each nucleotide and its relative frequency
            foreach (string seq in seqs)
            {
                
                string seqStr = seq;
                int count = 0;
                foreach (char c in seqStr.ToUpper(CI).ToCharArray()) // for a, then t, then c, then g, increment or add to that position
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
            } // end for each position. we have now constructed a dictionary


            // for each base position, sort and remove alleles with insufficient q
            for (int i = 0; i < freqList.Length; i++)
            {
                // Sort from most to least frequent
                freqList[i] = (from b in freqList[i] orderby b.Value descending select b)
                    .ToDictionary(pair => pair.Key, pair => pair.Value);

                double numBases = freqList[i].Values.Sum();
                char[] chars = freqList[i].Keys.ToArray();

                for (int k = 0; k < chars.Length; k++)
                {
                    // Convert to percentage and emove bases with insufficient qty for determination
                    char c = chars[k];
                    freqList[i][c] = (freqList[i][c] / numBases);

                    if (freqList[i][c] < 0.35) { freqList[i].Remove(c); }

                    if (freqList[i].Count == 0)
                    {
                        freqList[i].Add('N', 1);
                    }
                }

               
                foreach(KeyValuePair<char,double> p in freqList[i])
                {
                    if(alleles.ContainsKey(p.Key)){

                    }
                    if (p.Key == 'A' || p.Key == 'T' || p.Key == 'C' || p.Key == 'G')
                    {
                        if (masterFreqList[i] != null && masterFreqList[i].ContainsKey(p.Key))
                        {
                            masterFreqList[i][p.Key] += p.Value;
                        }
                        else
                        {
                            if (masterFreqList[i] == null)
                            {
                                masterFreqList[i] = new Dictionary<char, double>();
                            }
                            masterFreqList[i].Add(p.Key, p.Value);
                        }
                    }
                }
                
            }
            return freqList;
        }
        #endregion
        
        #endregion


        #region Private Static Methods

        #region helper methods

        /// <summary>
        /// Given a sequence, returns the reference ID, or null if sequences are unmapped
        /// </summary>
        /// <param name="sequences"></param>
        /// <returns></returns>
        private static string GetId(SAMAlignedSequence sequence)
        {
            if (!sequence.Flag.HasFlag(SAMFlags.UnmappedQuery))
            {
                return (sequence != null) ? sequence.RName : null;
            }
            else
            {
                return null;
            }
        }

        /// <summary>
        /// Given a SAMAlignedSequence, returns a string representation of the genetic sequence
        /// </summary>
        private static string GetSequence(SAMAlignedSequence seq)
        {
            String seqStr = seq.QuerySequence.ToString();
            return Regex.Split(seqStr, "\r\n")[0];
        }

        private static string GetReadQuality(SAMAlignedSequence seq)
        {
            String seqStr = seq.QuerySequence.ToString();
            return Regex.Split(seqStr, "\r\n")[1]; // todo aw something more efficient than regex?
        }

        // Get individual string
        private static string GetRgTag(SAMAlignedSequence seq)
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
        /// Add a sequence into a dictionary value item which represents a list of sequences. The list to which to add 
        /// to is found using key
        /// </summary>
        /// <param name="dict"></param>
        /// <param name="key"></param>
        /// <param name="seq"></param>
        private static void AddToDict(Dictionary<String, List<SAMAlignedSequence>> dict, string key, SAMAlignedSequence seq)
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
                throw new ArgumentException(Properties.Resources.INVALID_KEY);
            }
        }

        /// <summary>
        /// For a single sample in a single cluster, get top1, top2, ...topAll
        /// </summary>
        /// <param name="dict"></param>
        /// <returns></returns>
        private static List<double> GetSampleFrequencies(Dictionary<String, List<SAMAlignedSequence>> dict)
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

        /// <summary>
        /// Create a sequence dictionary from the current list of aligned sequences
        /// This is a dictionary where each distinct query sequence is represented as a key, and 
        /// a list of all individual sequences that share the same query string is stored
        /// as the value. The nested list in the dictionary is sorted in descending order of list size
        /// </summary>
        private static Dictionary<String, List<SAMAlignedSequence>> MakeSequenceDict(Collection<SAMAlignedSequence> seqs)
        {
            Dictionary<String, List<SAMAlignedSequence>> dict = new Dictionary<String, List<SAMAlignedSequence>>();

            foreach (SAMAlignedSequence seq in seqs)
            {
                AddToDict(dict, GetSequence(seq), seq);
            }
            return (from sequence in dict orderby sequence.Value.Count descending select sequence)
                    .ToDictionary(pair => pair.Key, pair => pair.Value);
        }

        #endregion

        /// <summary>
        /// Given an ordered collection of frequencies, return the sum of the top [ploidyLevel]
        /// </summary>
        /// <param name="ploidyLevel"></param>
        /// <param name="frequencies"></param>
        /// <returns></returns>
        private static double GetCountInPloidy(int ploidyLevel, Collection<double> frequencies)
        {
            double seqsCount = 0;
            int i = 0;
            foreach (double fx in frequencies)
            {
                if (i++ < ploidyLevel)
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

        #endregion


        #region Public Methods


        /// <summary>
        /// Calculate various metric values from the given list of sequences.
        /// </summary>
        /// <param name="clusterSequences">Sequences to add.</param>
        public void Calculate(Collection<SAMAlignedSequence> clusterSequences)
        {
            if (clusterSequences != null && clusterSequences.Count() > 0)
            {
                // Create various structures to store or index the sequence data in different ways
                this.sequencesThisCluster = clusterSequences;
                this.sequenceDict = MakeSequenceDict(new Collection<SAMAlignedSequence>(clusterSequences));
                MakeSampleDict();
                SetClusterReferenceIdAndSequence();

                // Get various counts once and once only, to prevent iterating through the above structures more than necessary
                countAll = clusterSequences.Count;
                countDistinct = sequenceDict.Count;
                countSamples = sampleDict.Count;
                populationPercentage = Math.Round(CountSamples / (double)numSamples, 2);

                // Create a more deeply nested structures
                MakeSampleSequenceDict();

                // Iterate through each structure the minimum number of times to calculate various things on the one loop
                // Various metrics are calculated from within these iterator methods
                IterateSequenceDict();  // alignment quality, read quality, frequencies of all distinct individuals
                IterateSampleDict();    // count how many sequences each person has in total (not distinct sequences)
                IterateSampleSequenceDict(); // cluster dirt, count distinct reads each individual has

                // Construct Phase input string for haplotyping
                ConstructPhaseSnpString();
            }
        }

        /// <summary>
        /// Return a string of tab-separated values for writing to file.
        /// </summary>
        /// <returns>A string of tab-separated values for writing to file.</returns>
        public override string ToString()
        {
            string header = (this.Id == "0") ? "#cluster_id\tcount_all_reads\tcount_distinct_reads\tnum_individuals\tdirt\talignment_qualities\tploidy_disagreement\tread_qualities\tpopulation_percentage\tavg_read_count_per_indiv_all\tavg_read_count_per_indiv_distinct\tnum_haplotypes" + Environment.NewLine : "";

            return header + 
                Id + "\t" + 
                CountAll + "\t" +
                CountDistinct + "\t" +
                CountSamples + "\t" +
                Dirt + "\t" +
                AlignmentQuality + "\t" +
                ploidyDisagreement + "\t" +
                ReadQuality + "\t" +
                PopulationPercentage + "\t" +
                Math.Round(SampleReadCountsAll.Average(), 2) + " \t " +
                Math.Round(SampleReadCountsDistinct.Average(), 2) + "\t" +
                numberOfHaplotypes;
            ;
        }

        /// <summary>
        /// Clear all structures
        /// </summary>
        public void Reset()
        {
            sequencesThisCluster.Clear();
            sequenceDict.Clear();
            sampleDict.Clear();
            sampleSequenceDict.Clear();
            readsInPloidyForIndividualsDict.Clear();
            readsNotInPloidyForIndividualsDict.Clear();
        }

        #endregion

    }
}
