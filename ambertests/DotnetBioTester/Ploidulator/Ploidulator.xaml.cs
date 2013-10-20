using Bio.Algorithms.Metric;
using Bio.IO.BAM;
using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Controls.DataVisualization;
using System.Windows.Controls.DataVisualization.Charting;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Animation;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Threading;

namespace Ploidulator
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {

        public delegate System.Delegate MyDelegate(string a, string b, string c, string d, string e, string f, bool? g);
        public delegate System.Delegate QuickDelegate();
        public delegate System.Delegate IntDelegate(int a);
        //public delegate System.Delegate SimpleDelegate(ref ClusterMetricHandlerPloidulator handler);
        ClusterMetricHandlerPloidulator handler = null;

        public MainWindow()
        {
            InitializeComponent();

            

        }

     

        // On click of GO button, begin ploidulation process
        private void Button_Click_Go(object sender, RoutedEventArgs e)
        {
            if (DataFileATextbox.Text != null)
            {
                // parse the user inputs here
                MyDelegate handler = ParseBAMMetric;
                handler.BeginInvoke(DataFileATextbox.Text, ExpectedPloidyTextbox.Text, 
                    DirtCutoffTextbox.Text, AlignmentQualCutoffTextbox.Text,
                    ReadQualCutoffTextbox.Text, NumSamplesTextbox.Text, OutputToFile.IsChecked, null, null);
            }
        }


        // Demo that window thread is live
        /*private void Button_Click_Circle(object sender, RoutedEventArgs e)
        {
            Dispatcher.BeginInvoke(
                System.Windows.Threading.DispatcherPriority.Normal,
                new OneArgDelegate(DrawCircle),
                "");
        }*/

        // Open file selector window and fetch selected file path
        private void Button_Click_GetFile(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog dialog = new Microsoft.Win32.OpenFileDialog();
            dialog.DefaultExt = ".bam";
            dialog.Filter = "BAM files (*.bam)|*.bam";
            if ((bool)dialog.ShowDialog())
            {
                DataFileATextbox.Text = dialog.FileName;
            } 
        }

        private System.Delegate UpdateGui()
        {
            DataFileATextbox.Text = "";
            FileExpander.IsExpanded = false;
            return null;
        }

        private System.Delegate SetProgressBar(int max)
        {
            ProgressBar.IsEnabled = true;
                
            ProgressBar.Minimum = 0;

            ProgressBar.Maximum = max;
            ProgressBar.IsIndeterminate = false;
            ProgressBar.ToolTip = "Found " + max.ToString() + "clusters to analyse ...";
            ProgressBar.Visibility = System.Windows.Visibility.Visible;
                
            return null;
        }

        // Draw a simple non-functional circle for testing
        public System.Delegate DrawCircle(string message)
        {
            Ellipse circle = new Ellipse();
            SolidColorBrush mySolidColorBrush = new SolidColorBrush();
            mySolidColorBrush.Color = (Color)ColorConverter.ConvertFromString("#FFFFFF00");
            circle.Fill = mySolidColorBrush;
            circle.StrokeThickness = 2;
            circle.Stroke = Brushes.Black;
            circle.Width = 75;
            circle.Height = 75;
            wpMain.Children.Add(circle);
            return null;
        }
        /*private System.Delegate LoadMorePies(ref ClusterMetricHandlerPloidulator handler)
        {
            if(handler != null)
            {
                handler.MaxNumCharts = handler.MaxNumCharts + 10;
            }
            
            return null;
        }*/

        private System.Delegate ParseBAMMetric(string filename, string ploidy, string dirtCutoff, string alignQualCutoff, 
            string readQualCutoff, string numSamples, bool? outputToFile)
        {
            // Check user input
            string newName = filename.Split(new char[]{'.'})[0];

            // Set up the handler
            handler = new ClusterMetricHandlerPloidulator(newName + "_pl", Convert.ToInt32(ploidy), Convert.ToDouble(dirtCutoff),
                Convert.ToInt32(alignQualCutoff), readQualCutoff, Convert.ToInt32(numSamples), outputToFile, Dispatcher, wpMain); 
            //handler.WriteToFile(false); // todo oops, it still opens file
            
            // Start the parser
            BAMParser parser = new BAMParser();
            int clusters = 0;
            using (Stream readStream = new FileStream(filename, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                SAMAlignmentHeader header = parser.GetHeader(readStream);
                handler.Header = header;
                clusters = header.ReferenceSequences.Count();
                Console.WriteLine("total num seq" + clusters);
            }

            Dispatcher.BeginInvoke(
                System.Windows.Threading.DispatcherPriority.Normal,
                new IntDelegate(SetProgressBar), clusters);

            // Everything appears ok, so close the GUI window
            Dispatcher.BeginInvoke(
                System.Windows.Threading.DispatcherPriority.Normal,
                new QuickDelegate(UpdateGui));

            // Pass each sequence to the handler to be processed
            

            double increment = 1; // increase this if we only want to update the progress bar every x many clusters
            int updFor = -1;
            int count = -1;
            foreach (SAMAlignedSequence se in parser.ParseSequence(filename))
            {
                handler.Add(se);
                count = handler.ClusterCount;
                if ((double)count % increment == 0 && count != updFor)
                {
                    updFor = count;
                    Console.WriteLine("DIVISIBLE");
                    Dispatcher.BeginInvoke(
                        System.Windows.Threading.DispatcherPriority.Normal,
                        new IntDelegate(UpdateProgress), handler.ClusterCount);
                }
                
                
            }

            // Close handler and return
            handler.SetComplete();
            handler.Dispose();
            parser.Dispose();
            
            return null;
        }

        private System.Delegate UpdateProgress(int index)
        {
            //index = index * 100;

            if(index >= ProgressBar.Minimum && index <= ProgressBar.Maximum)
            {
                ProgressBar.Value = index;
                ProgressBar.ToolTip = "Analysed " + index.ToString() + " out of " + ProgressBar.Maximum.ToString() + " clusters";
            }
            else
            {
                throw new Exception("Invalid progress bar value");
            }
            if (ProgressBar.Maximum == index)
            {
                ProgressBar.Visibility = System.Windows.Visibility.Visible;
            }
            return null;
            
        }

        private void ProgressBar_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {

        }
        

        /*private void Button_Click_LoadMore(object sender, RoutedEventArgs e)
        {
            if (this.handler != null)
            {
                SimpleDelegate del = LoadMorePies;
                del.BeginInvoke(ref handler, null, null);
            }    
        }*/


    }
}
