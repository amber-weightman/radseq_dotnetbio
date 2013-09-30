using Bio;
using Bio.Algorithms.Metric;
using Bio.IO.BAM;
using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DotnetBioTester
{
    class Program
    {
        static void Main(string[] args)
        {
            string fileName = @"E:\Harvard\final-all_bam-merged.bam";
            ParseBAMMetric(fileName);
            //Dirt();

            //int numClusters = header.ReferenceSequences.Count;

        }


        static void Dirt()
        {
               //load the header for the BAM File
             //string fname = @"C:\Users\Nigel\SkyDrive\Ploidy\final-all_bam-merged.bam" ;
            string fname = @"E:\Harvard\final-all_bam-merged.bam";
             //fname = "../../../../../../../../Harvard/final-all_bam-merged.bam";

                      BAMParser parser = new BAMParser();
                      var temp=parser.ParseRange (fname, 0);
                      var header=temp.Header;
                      int numClusters = header.ReferenceSequences.Count;
                      //now run through each cluster and calculate cdirt in parrallel
                      var cdirts= Enumerable.Range (0, numClusters).AsParallel().AsOrdered().Select (refID =>
                     {
                                     
                            BAMParser contigParser = new BAMParser ();
                            var samSeqs = contigParser.ParseRange (fname, refID);
                            //Convert sequences in to smaller types and group them by sample
                            var seqs = samSeqs.QuerySequences.Select (x =>
                                                new { Sample=x.OptionalFields.Where (z => z.Tag == "RG" ).Select (y => y.Value).First (),
                                                      Sequence=(x.QuerySequence as QualitativeSequence ).ToString ().Substring (0, (int)x.QuerySequence.Count)});
                            var seqsBySample = seqs.GroupBy(x => x.Sample);
                            //calculate the number of reads in the top two and the remaining two for each sample
                            List< Tuple< double, double>> results = new List<Tuple <double , double >> ();
                            foreach ( var sample in seqsBySample) {
                                   var counts = sample.GroupBy (z => z.Sequence).Select (x => x.Count ()).ToList ();
                                  counts.Sort ((x,y) => -x.CompareTo (y));
                                  results.Add ( new Tuple< double, double> ((double)counts.Take (2).Sum (),
                                                                           (double)counts.Skip (2).Sum ()));
                           }
                            //now calculate the c-dirt
                            double totalReads = results.Select (x => x.Item1 + x.Item2).Sum ();
                            double minorityReads = results.Select (x => x.Item2).Sum();
                            double cdirt = minorityReads / totalReads;
                            return new Tuple< int, double> (refID, cdirt);
                     }).ToList();
                      //now just output
                      
            string writeF = "../../../../../../../../Harvard/final-all_bam-merged.csv";
              //StreamWriter sw = new StreamWriter(@"C:\Users\Nigel\SkyDrive\Ploidy\CdirtND.csv" );
              StreamWriter sw = new StreamWriter(writeF);
            //sw.WriteLine(x.Item1.ToString()+"," +x.Item2.ToString());

            cdirts.ForEach(x => Console.WriteLine(x.Item1.ToString() + "," + x.Item2.ToString()));
                     cdirts.ForEach(x=> sw.WriteLine(x.Item1.ToString()+"," +x.Item2.ToString()));
                     sw.Close ();
                                  
              }
        

    
        static void ParseBAMMetric(string filename)
        {
            // Where a different metric handler would perform different operations on each alignment/cluster
            ClusterMetricHandler handler = new ClusterMetricHandler("../../../../../../../../Harvard/final-all_bam-merged");

            // False also applies by default. If set to true, reads will be written to a new file
            // (sam or bam), unless they are filtered out by the handler (i.e. cluster dirt outside range)
            handler.WriteToFile(false); // todo oops, it still opens file

            // Where "false" indicates we want the parser to execute the metric operation only
            // -- we do not want to read BAM file into memory
            // ("true" would make the parser do the same thing it has always done (overriding "Lite"),
            // with the metric operation being performed in addition. Omitting the bool would utilise
            // default memory storage, which currently stores the SequenceAlignmentMap in memory (same as passing "true"))
            BAMParser parser = new BAMParser(handler, false);
            Console.WriteLine("parser created, about to parse file and create .metr metric...");

            // It still returns a sequence alignment map, but in this case the map will be null
            SequenceAlignmentMap pp = parser.Parse(filename);
        }

        static SequenceAlignmentMap ParseBAM(string filename)
        {
            // Where the parser does the same thing it has always done, even though I refactored it.
            // The existing tests still need to be run to make sure I didn't break anything expected
            BAMParser parser = new BAMParser();
            Console.WriteLine("parser created, about to parse file into memory...");
            SequenceAlignmentMap pp = parser.Parse(filename);
            return pp;
        }


    }


}
