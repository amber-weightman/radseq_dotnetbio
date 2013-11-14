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

        private static CultureInfo ci = new CultureInfo("en-AU");

        #endregion

        #region Private Fields

        /// <summary>
        /// The expected ploidy level of the organism under study (2 by default)
        /// </summary>
        private int expectedPloidy = 2;

        /// <summary>
        /// The number of samples represented in the entire population
        /// </summary>
        private int numSamples;

        /// <summary>
        /// The number of haplotypes in this cluster
        /// </summary>
        private int numberOfHaplotypes = -1;

        #region storing sequences

        /// <summary>
        /// All reads which are in-ploidy for the owning individual
        /// </summary>
        private Collection<SAMAlignedSequence> readsInPloidyForIndividuals = null;

        /// <summary>
        /// All reads which are not in-ploidy for the owning individual
        /// </summary>
        private Collection<SAMAlignedSequence> readsNotInPloidyForIndividuals = null;

        /// <summary>
        /// All reads which are in-ploidy for the owning individual, as a dictionary. 
        /// The nested list in the dictionary is sorted in descending order of list size
        /// </summary>
        private Dictionary<string, List<SAMAlignedSequence>> readsInPloidyForIndividualsDict = null;

        /// <summary>
        /// All reads which are not in-ploidy for the owning individual, as a dictionary.
        /// The nested list in the dictionary is sorted in descending order of list size
        /// </summary>
        private Dictionary<string, List<SAMAlignedSequence>> readsNotInPloidyForIndividualsDict = null;

        /// <summary>
        /// A dictionary where each distinct query sequence is represented as a key, and 
        /// a list of all individual sequences that share the same query string is stored
        /// as the value
        /// The nested list in the dictionary is sorted in descending order of list size
        /// </summary>
        private Dictionary<String, List<SAMAlignedSequence>> sequenceDict = null;

        /// <summary>
        /// A dictionary where each sample individual is represented as a key, and 
        /// a list of all individual sequences for that individual is stored
        /// as the value
        /// The nested list in the dictionary is sorted in descending order of list size
        /// </summary>
        private Dictionary<String, List<SAMAlignedSequence>> sampleDict = null;

        /// <summary>
        /// Experimental
        /// </summary>
        private Dictionary<String, Dictionary<String, List<SAMAlignedSequence>>> sampleSequenceDict = null;

        /// <summary>
        /// List of all sequences.
        /// </summary>
        private Collection<SAMAlignedSequence> sequences = null;

        #endregion

        /// <summary>
        /// Cluster ID (unique ID of reference sequence against which all sequences in the cluster
        /// are aligned)
        /// </summary>
        private string id;

        /// <summary>
        /// Reference sequence for the cluster (sequence which all reads are aligned to)
        /// </summary>
        private string referenceSequence;

        /// <summary>
        /// Formatted string as required as input by PHASE. For each allele position, 'S' indicates biallelic and
        /// 'M' indicates multiallelic 
        /// </summary>
        private string loci;

        /// <summary>
        /// Formatted string as required as input by PHASE. Multi-line string where probable genotype for each
        /// sequence appears on two lines per sequence
        /// </summary>
        private string phaseData;

        /// <summary>
        /// Total number of all sequences in the cluster
        /// </summary>
        private int countAll;

        /// <summary>
        /// Total number of distinct sequences in the cluster
        /// </summary>
        private int countDistinct;

        /// <summary>
        /// Total number of samples (individuals) represented in the cluster
        /// </summary>
        private int countSamples;

        #endregion

        #region Constructors

        /// <summary>
        /// The default constructor.
        /// </summary>
        public ClusterMetric()
        {
        }

        /// <summary>
        /// Non-default constructor
        /// </summary>
        /// <param name="ploidy">Expected ploidy level of the organism</param>
        /// <param name="samples">Number of samples in data file</param>
        public ClusterMetric(int ploidy, int samples)
        {
            expectedPloidy = ploidy;
            numSamples = samples;
        }

        #endregion

        #region Properties

        /// <summary>
        /// Cluster ID (unique ID of reference sequence against which all sequences in the cluster
        /// are aligned)
        /// </summary>
        public string Id { get { return id; } }

        /// <summary>
        /// Formatted string as required as input by PHASE. For each allele position, 'S' indicates biallelic and
        /// 'M' indicates multiallelic 
        /// </summary>
        public string PhaseLoci { get { return loci; } }

        /// <summary>
        /// Formatted string as required as input by PHASE. Multi-line string where probable genotype for each
        /// sequence appears on two lines per sequence
        /// </summary>
        public string PhaseData { get { return phaseData; } }

        /// <summary>
        /// The number of probable haplotypes in this cluster
        /// </summary>
        public int NumberOfHaplotypes { get { return numberOfHaplotypes; } set { numberOfHaplotypes = value; } }

        /// <summary>
        /// Total number of all sequences in the cluster
        /// </summary>
        public int CountAll { get { return countAll; } }

        

        /// <summary>
        /// Total number of distinct sequences in the cluster
        /// </summary>
        public int CountDistinct { get { return countDistinct; } }

        /// <summary>
        /// Total number of samples (individuals) represented in the cluster
        /// </summary>
        public int CountSamples { get { return countSamples; } }


        /// <summary>
        /// A dictionary where each distinct query sequence is represented as a key, and 
        /// a list of all individual sequences that share the same query string is stored
        /// as the value. 
        /// </summary>
        [System.Diagnostics.CodeAnalysis.SuppressMessage("Microsoft.Design", "CA1006:DoNotNestGenericTypesInMemberSignatures")]
        public Dictionary<String, List<SAMAlignedSequence>> SequenceDictionary { get { return sequenceDict; } }

        /// <summary>
        /// A dictionary where each sample individual is represented as a key, and 
        /// a list of all individual sequences for that individual is stored
        /// as the value
        /// </summary>
        [System.Diagnostics.CodeAnalysis.SuppressMessage("Microsoft.Design", "CA1006:DoNotNestGenericTypesInMemberSignatures")]
        public Dictionary<String, List<SAMAlignedSequence>> SampleDictionary { get { return sampleDict; } }

        /// <summary>
        /// List of all sequences.
        /// </summary>
        public Collection<SAMAlignedSequence> Sequences { get { return sequences; } }



        /// <summary>
        /// Frequency distribution of all distinct reads, regardless of individual
        /// </summary>
        private Collection<int> frequencyDistributionSequences;

      
        /// <summary>
        /// Average frequencies for each sequence, calculated per sample per cluster and averaged
        /// out to cluster
        /// </summary>
        public Collection<double> ClusterSequenceFrequencies { get { return individualSequenceDistributions; } }

        /// <summary>
        /// Some measure of alignment quality using MAPQ
        /// </summary>
        public Collection<double> ClusterAlignmentQualities { get { return clusterAlignmentQualities; } }

        /// <summary>
        /// Placeholder (determine if cluster is good based on metrics etc)
        /// if the metric determines that we want to keep this read for downstream analysis
        /// </summary>
        public bool Good { get { return good; } set { this.good = value; } }

        
        private bool good = false;

        public double Dirt { get { return dirt; } }
        private double dirt;

        /// <summary>
        /// being the avg qual for distinct reads per indiv which are within ploidy (the reads on the right side of cluster dirt)
        /// </summary>
        public double AlignmentQuality { get { return alignmentQuality; } }
        private double alignmentQuality;

        /// <summary>
        /// Are all individuals represented in the cluster? This is the percentage of whole-population representation
        /// </summary>
        public double PopulationPercentage { get { return populationPercentage; } }
        private double populationPercentage;

        public Collection<double> SampleReadCountsAll { get { return sampleReadCountsAll; } }
        private Collection<double> sampleReadCountsAll;

        public Collection<double> SampleReadCountsDistinct { get { return sampleReadCountsDistinct; } }
        private Collection<double> sampleReadCountsDistinct;


        public double ReadQuality { get { return readQuality; } }
        private double readQuality;

        #endregion

        #region Public Methods

        /// <summary>
        /// String of tab-separated values for writing to file by MetricFormatter.
        /// todo fixme we want the MetricFormater to decide how to write each line, but
        /// for now until I know what values are being written I will let the Metric handle this.
        /// </summary>
        public override string ToString()
        {
            string header = (this.Id == "0") ? "#cluster_id\tcount_all_reads\tcount_distinct_reads\tnum_individuals\tdirt\talignment_qualities_all\tploidy_disagreement_unnormalised\tread_qualities_all\tpopulation_percentage\tavg_read_count_per_indiv_all\tavg_read_count_per_indiv_distinct\tnum_haplotypes" + Environment.NewLine : "";

            return header  + Id + "\t" + CountAll + "\t" +
                CountDistinct + "\t" +
                CountSamples + "\t" + 
                Dirt + "\t" + 
                AlignmentQuality  + "\t(" + 
                ploidyDisagreement + ")\t" +
                ReadQuality + "\t" + 
                PopulationPercentage + "\t" + 
                Math.Round(SampleReadCountsAll.Average(), 2) + " \t " + 
                Math.Round(SampleReadCountsDistinct.Average(), 2)  + "\t" + 
                numberOfHaplotypes;
            ;
        }
        /// <summary>
        /// List of frequencies for each cluster, from most to least frequent (sums to 1). Ideally top [expectedPloidy] sequences
        /// would together be ~1 (100%)
        /// </summary>
        private Collection<double> individualSequenceDistributions = null;

        /// <summary>
        /// The average alignment quality of each (distinct??) read in this cluster
        /// </summary>
        private Collection<double> clusterAlignmentQualities = null;



        /// <summary>
        /// Calculate metric values from the given list of sequences.
        /// </summary>
        /// <param name="clusterSequences">Sequences to add.</param>

        public void Calculate(Collection<SAMAlignedSequence> clusterSequences)
        {
            // Create various structures to store or index the sequence data in different ways
            this.sequences = clusterSequences;
            sequenceDict = MakeSequenceDict(new Collection<SAMAlignedSequence>(clusterSequences));
            sampleDict = MakeSampleDict(clusterSequences);
            sampleSequenceDict = MakeNestedSequenceDict(sampleDict);
            
            // Get the counts once and once only
            countAll = clusterSequences.Count;
            countDistinct = sequenceDict.Count;
            countSamples = sampleDict.Count;
            populationPercentage = Math.Round(CountSamples / (double)numSamples, 2);

            SetClusterReferenceIdAndSequence(); // has nothing to do with iterating
            ConstructPhaseMasterTemplate();
            
            // Iterate through each structure the minimum number of times to calculate various things on the one loop
            IterateSequenceDict(); // does alignment and read quality and frequencies of all distinct individuals
            IterateSampleDict();    // gets a count of how many sequences each person has in total (not distinct sequences)
            IterateSampleSequenceDict();
            
            ploidyDisagreement = FindPloidyDisagreement();
            
            // At this point every value should be set and nothing more should need to be done
        }

        #region iterateSequences

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


        private void SetClusterReferenceIdAndSequence()
        {
            id = GetId(sequences[0]); 
            referenceSequence = GetSequence(sequences[0]);
        }


       

        private void ConstructPhaseMasterTemplate()
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

        // precondition seqs is an ordered list
        private Dictionary<char, double>[] BaseFrequencies(string[] seqs, ref Dictionary<char, double>[] masterFreqList)
        {
            
            // A dictionary of char:double for each base pair position
            Dictionary<char, double>[] freqList = new Dictionary<char, double>[seqs[0].Length];

            //For each position, get each nucleotide and its relative frequency
            foreach (string seq in seqs)
            {
                
                string seqStr = seq;
                int count = 0;
                foreach (char c in seqStr.ToUpper(ci).ToCharArray()) // for a, then t, then c, then g, increment or add to that position
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

        #region anySequenceList


        #endregion




        #region IterateSequenceDict

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

            frequencyDistributionSequences = frequencies;
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

     

     

        private double FindPloidyDisagreement()
        {
            return Math.Round(((readsInPloidyForIndividualsDict.Count 
                + readsNotInPloidyForIndividualsDict.Count) - CountDistinct) / (double)CountDistinct, 2);
        }
        private double ploidyDisagreement;


        #endregion

        
       /// <summary>
       /// Sets, among other things, dirt. Also produces a count of the number of distinct reads each individual has
       /// </summary>
        private void IterateSampleSequenceDict()
        {
            individualSequenceDistributions = SampleFrequenciesAvg();
            dirt = Math.Round(1 - GetCountInPloidy(expectedPloidy, individualSequenceDistributions), 2);
            // dirt = readsNotInPloidyForIndividualsDict.Count / (double)(readsInPloidyForIndividualsDict.Count + readsNotInPloidyForIndividualsDict.Count);

            sampleReadCountsDistinct = new Collection<double>();
            foreach (Dictionary<string, List<SAMAlignedSequence>> seqList in sampleSequenceDict.Values)
            {
                sampleReadCountsDistinct.Add(seqList.Count);
            }   
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
            readsInPloidyForIndividuals.Clear();
            readsNotInPloidyForIndividuals.Clear();

            readsInPloidyForIndividualsDict.Clear();
            readsNotInPloidyForIndividualsDict.Clear();
        }

        #endregion


        #region Private Static Methods



        


        private static double GetCountInPloidy(int ploidyLevel, Collection<double> frequencies)
        {
            int count = 0;
            double seqsCount = 0;
            foreach (double fx in frequencies)
            {
                if (count++ < ploidyLevel)
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

        /// <summary>
        /// Create a sample dictionary from the current list of aligned sequences
        /// This is a dictionary where each sample individual is represented as a key, and 
        /// a list of all individual sequences for that individual is stored
        /// as the value
        /// Each value in the dictionary is a shallow copy of an aligned sequence in sequences
        /// </summary>
        private static Dictionary<String, List<SAMAlignedSequence>> MakeSampleDict(Collection<SAMAlignedSequence> seqs)
        {
            Dictionary<String, List<SAMAlignedSequence>> dict = new Dictionary<String, List<SAMAlignedSequence>>();

            foreach (SAMAlignedSequence seq in seqs)
            {
                AddToDict(dict, GetRgTag(seq), seq);
            }
            return dict;
        }

        // Add sequence to dictionary using key
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
                throw new ArgumentException("Invalid dictionary key");
            }
        }


        private Collection<double> SampleFrequenciesAvg()
        {

            List<double>[] thisClusterSampleFrequencies = new List<double>[sampleSequenceDict.Count];
            int i = 0;
            int maxlen = 0;
            foreach (Dictionary<String, List<SAMAlignedSequence>> sample in sampleSequenceDict.Values) // for each indiv
            {
                thisClusterSampleFrequencies[i] = GetSampleFrequencies(sample);
                maxlen = (thisClusterSampleFrequencies[i].Count > maxlen) ? thisClusterSampleFrequencies[i].Count : maxlen;
                i++;
            }

            //int maxlen = 0;
            /*foreach (List<double> fx in thisClusterSampleFrequencies)
            {
                maxlen = (fx.Count > maxlen) ? fx.Count : maxlen;
            }*/

            double[] temp = new double[maxlen];
            Array.Clear(temp, 0, maxlen); // fill with 0s
            foreach (List<double> fx in thisClusterSampleFrequencies)
            {
                int j;
                for (j = 0; j < fx.Count; j++)
                {
                    temp[j] += fx[j];
                }
            }
            for (int j = 0; j < temp.Length; j++)
            {
                temp[j] = temp[j] / (double)thisClusterSampleFrequencies.Length; // divide by num samples
            }

            Collection<double> d = new Collection<double>();
            foreach (double t in temp)
            {
                d.Add(t);
            }

            return d;
        }



        /// <summary>
        /// For a single sample in a single cluster, get top1, top2, ...topAll
        /// </summary>
        /// <param name="dict"></param>
        /// <returns></returns>
        private  List<double> GetSampleFrequencies(Dictionary<String, List<SAMAlignedSequence>> dict)
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


        #region Private Methods

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
                superDict[entry.Key] = MakeSequenceDict(new Collection<SAMAlignedSequence>(entry.Value)); // superDict["sampleName"] = [a sequence dictionary]
                if (readsInPloidyForIndividuals == null)
                {
                    readsInPloidyForIndividuals = new Collection<SAMAlignedSequence>();
                    readsNotInPloidyForIndividuals = new Collection<SAMAlignedSequence>();
                }

                int i = 0;
                foreach (List<SAMAlignedSequence> t in superDict[entry.Key].Values) // for each sequence dictionary
                {
                    if(i++ < expectedPloidy){
                        foreach(SAMAlignedSequence s in t)
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
            Debug.Assert((readsInPloidyForIndividuals.Count + readsNotInPloidyForIndividuals.Count) == sequences.Count);
            

            // todo fixme there might be a sequence in both buckets that shares the same query string

            return superDict;
        }

        #endregion

    }
    

}
