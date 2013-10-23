using Bio.Algorithms.Metric;
using Bio.IO.BAM;
using Bio.IO.SAM;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
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

        public delegate System.Delegate MyDelegate(string a, string b, string c, string d, string e, string f, string f2,
            bool? g, bool? h, bool? i, bool? j, bool? k);
        public delegate System.Delegate QuickDelegate();
        public delegate System.Delegate IntDelegate(int a);
        public delegate System.Delegate StatsDelegate(int a, int b, int c, double d, double e, double f, double g, double h, double i, double j, double k);
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
                    ReadQualCutoffTextbox.Text, PopPercentTextbox.Text, NumSamplesTextbox.Text, OutputToFile.IsChecked, 
                    MetricFileParent.IsChecked, MetricFileChild.IsChecked, OutputOverviewParent.IsChecked, OutputOverviewChild.IsChecked,
                    null, null);
            }
        }

        private void Button_Click_Abort(object sender, RoutedEventArgs e)
        {
            AbortingMessage.Visibility = System.Windows.Visibility.Visible;
            handler.Abort();
        }

        private void Button_Click_Preview(object sender, RoutedEventArgs e)
        {
            // not implemented
        }



        // Open file selector window and fetch selected file path
        private void Button_Click_GetFile(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog dialog = new Microsoft.Win32.OpenFileDialog();
            dialog.DefaultExt = ".bam";
            dialog.Filter = "BAM files (*.bam)|*.bam";
            if ((bool)dialog.ShowDialog())
            {
                DataFileATextbox.Text = dialog.FileName;
                TimerLabel.Visibility = System.Windows.Visibility.Hidden;
            } 
        }
        private System.Delegate ProcessInputSequences()
        {
            if (inputQueue.Count == 0)
            {
                Thread.Sleep(5000); // sleep 5 seconds
            }
            else if (inputQueue.Count > 100)
            {
                Thread.Sleep(20000); // sleep 20 seconds
            }
            else
            {
                // do processing

            }
            return null;
        }

        private System.Delegate UpdateGui_BeganParsing()
        {
            isProcessingFile = true;

            // Disable running new job
            ToggleSearchable(false);

            // Show stats panel
            StatsPanel.Visibility = System.Windows.Visibility.Visible;

            // Show progress bar
            ProgressStatusBar.Visibility = System.Windows.Visibility.Visible;
            AbortButton.Visibility = System.Windows.Visibility.Visible;

            MenuExpander0.IsExpanded = true;
            MenuExpander1.IsExpanded = true;
            MenuExpander2.IsExpanded = false;
            

            StartTimer();

            return null;
        }

        private System.Windows.Threading.DispatcherTimer dispatcherTimer;
        private void StartTimer()
        {
            TimerLabel.Content = "";
            dispatcherTimer = new System.Windows.Threading.DispatcherTimer();
            dispatcherTimer.Tick += new EventHandler(dispatcherTimer_Tick);
            dispatcherTimer.Interval = new TimeSpan(0, 0, 1);
            dispatcherTimer.Start();
            startedTime = DateTime.Now;
            TimerLabel.Visibility = System.Windows.Visibility.Visible;
        }

        private void StopTimer()
        {
            dispatcherTimer.Stop();
            //TimerLabel.Visibility = System.Windows.Visibility.Hidden;
        }
        private DateTime startedTime;

        private void dispatcherTimer_Tick(object sender, EventArgs e)
        {
            TimeSpan time = DateTime.Now - startedTime;
            TimerLabel.Content = "Time elapsed: " + time.Hours.ToString() + ":" + time.Minutes.ToString("00") + ":" + time.Seconds.ToString("00");
            // todo probably a better way?
        }

        private bool isProcessingFile = false;
        private System.Delegate UpdateGui_FinishedParsing()
        {
            // Enable running new job
            ToggleSearchable(true);

            // Hide progress bar
            ProgressStatusBar.Visibility = System.Windows.Visibility.Hidden;
            AbortButton.Visibility = System.Windows.Visibility.Hidden;

            // Hide any messages
            AbortingMessage.Visibility = System.Windows.Visibility.Hidden;

            // Reset to black
            SamplesDisplay.Foreground = Brushes.Black;
            SamplesDisplayLabel.Foreground = Brushes.Black;

            StopTimer();
            isProcessingFile = false;
            return null;
        }

        private void ToggleSearchable(bool isSearchable)
        {
            FileExpander.IsExpanded = isSearchable;
            OutputFileExpander.IsExpanded = false;
            DataFileAButton.IsEnabled = isSearchable;
            DataFileATextbox.IsEnabled = isSearchable;
            ExpectedPloidyTextbox.IsEnabled = isSearchable;
            AlignmentQualCutoffTextbox.IsEnabled = isSearchable;
            DirtCutoffTextbox.IsEnabled = isSearchable;
            NumSamplesTextbox.IsEnabled = isSearchable;
            ReadQualCutoffTextbox.IsEnabled = isSearchable;
            PopPercentTextbox.IsEnabled = isSearchable;
            LaunchButton.IsEnabled = isSearchable;
        }

        
        private System.Delegate UpdateStatsPanel(int numGoodClusters, int numClustersParsed, int maxSampleCount,
            double maxMapq, double maxReadq, double avgDirt, double avgMapq, double avgReadq, 
            double avgDirtGood, double avgMapqGood, double avgReadqGood)
        {
            if (numClustersParsed != 0)
            {


                UpdateProgress(numClustersParsed);

                Console.Write("-begin upd stats panel-");
                SamplesDisplay.Content = maxSampleCount.ToString();
                MAPQdisplay.Content = maxMapq.ToString();
                READQdisplay.Content = maxReadq.ToString();

                DirtDisplay.Content = avgDirt.ToString();
                AvgMAPQDisplay.Content = avgMapq.ToString();
                AvgREADQDisplay.Content = avgReadq.ToString();

                DirtDisplayGood.Content = avgDirtGood.ToString();
                AvgMAPQDisplayGood.Content = avgMapqGood.ToString();
                AvgREADQDisplayGood.Content = avgReadqGood.ToString();


                double numGood = Math.Round(numGoodClusters / (double)numClustersParsed * 100, 2);
                double numBad = 100 - numGood;
                GoodBadProgressBar.Value = numGood;
                GoodCountLabel.Content = numGood + "% good";
                BadCountLabel.Content = numBad + "% bad";

                string plural = numGoodClusters == 1 ? "" : "s";
                string plural2 = (numClustersParsed - numGoodClusters) == 1 ? "" : "s";

                GoodCountLabel.ToolTip = numGood + "% of clusters are good (" + numGoodClusters + " cluster" + plural + ")";
                BadCountLabel.ToolTip = numBad + "% of clusters are bad (" + (numClustersParsed - numGoodClusters) + " cluster" + plural2 + ")";
                Console.Write("-end upd panel-");

                int a = Convert.ToInt32(maxSampleCount);
                int b = Convert.ToInt32(NumSamplesTextbox.Text);
                if ((int)a < (int)b && (int)a != 0)
                {
                    Console.Write("WRONG NUMBER OF SAMPLES");
                    SamplesDisplay.Foreground = Brushes.Red;
                    SamplesDisplayLabel.Foreground = Brushes.Red;
                }
                else
                {
                    SamplesDisplay.Foreground = Brushes.Black;
                    SamplesDisplayLabel.Foreground = Brushes.Black;
                    Console.Write("num samples ok for " + a + " and " + b);

                }

                if (numGoodClusters > 0)
                {
                    MenuExpander2.IsExpanded = true;
                }

            }
            return null;
        }

        private System.Delegate SetProgressBar(int max)
        {
            ProgressBar.IsEnabled = true;
            ProgressBar.Minimum = 0;
            ProgressBar.Maximum = max;
            ProgressBar.IsIndeterminate = false;
            ProgressBar.ToolTip = "Found " + max.ToString() + " clusters to analyse ...";
            ProgressBar.Visibility = System.Windows.Visibility.Visible;
            AbortButton.Visibility = System.Windows.Visibility.Visible;
            LoadingBarLabel.Visibility = System.Windows.Visibility.Visible;
                
            return null;
        }


        private System.Delegate ParseBAMMetric(string filename, string ploidy, string dirtCutoff, string alignQualCutoff, 
            string readQualCutoff, string popPercent, string numSamples, bool? outputToFile, bool? metricFileParent, 
            bool? metricFileChild, bool? outputOverviewParent, bool? outputOverviewChild)
        {
            // Check user input
            string newName = filename.Split(new char[]{'.'})[0];

            // Set up the handler
            handler = new ClusterMetricHandlerPloidulator(newName + "_pl", Convert.ToInt32(ploidy), Convert.ToDouble(dirtCutoff),
                Convert.ToDouble(alignQualCutoff), Convert.ToDouble(readQualCutoff), Convert.ToDouble(popPercent), Convert.ToInt32(numSamples), outputToFile, Dispatcher, wpMain);

            handler.WriteToFilteredBam = outputToFile == true;
            handler.WriteClusterMetricOriginal = metricFileParent == true;
            handler.WriteClusterMetricFiltered = metricFileChild == true;
            handler.WriteOverviewMetricOriginal = outputOverviewParent == true;
            handler.WriteOverviewMetricFiltered = outputOverviewChild == true;


            // Start the parser
            BAMParser parser = new BAMParser();
            int clusters = 0;
            using (Stream readStream = new FileStream(filename, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                SAMAlignmentHeader header = parser.GetHeader(readStream);
                handler.Header = header;
                clusters = header.ReferenceSequences.Count();
            }

            Dispatcher.BeginInvoke(
                System.Windows.Threading.DispatcherPriority.Normal,
                new IntDelegate(SetProgressBar), clusters);

            // Everything appears ok, so close the GUI window
            Dispatcher.BeginInvoke(
                System.Windows.Threading.DispatcherPriority.Normal,
                new QuickDelegate(UpdateGui_BeganParsing));



            // Pass each sequence to the handler to be processed
            double incrementOne = 1; // the progress bar will update every [incrementOne] many clusters
            double incrementTwo = 5; // the progress bar will update every [incrementTwo] many clusters
            int increaseIncrementAfter = 100; // after this many clusters, increase the increment at which we update the gui
            double increment;

            int updFor = -1, count = -1; // whether we have already updated the gui for this cluster

            inputQueue = new Queue<SAMAlignedSequence>();

            // new thread which sleeps and checks whether inpuq q has items
            Dispatcher.BeginInvoke(
                System.Windows.Threading.DispatcherPriority.Normal,
                new QuickDelegate(ProcessInputSequences));

            /*foreach (SAMAlignedSequence se in parser.ParseSequence(filename))
            {
                inputQueue.Enqueue(se);
            }*/

                foreach (SAMAlignedSequence se in parser.ParseSequence(filename))
                {
                    if (handler.Add(se))
                    {
                        count = handler.ClusterCount; 
                        increment = (count < increaseIncrementAfter) ? incrementOne : incrementTwo;
                        if ((count == 1 || (double)count % increment == 0) && count != updFor)
                        {
                            updFor = count;

                            Dispatcher.BeginInvoke(
                                System.Windows.Threading.DispatcherPriority.Normal,
                                new StatsDelegate(UpdateStatsPanel), handler.GoodCount, count,
                                    handler.MaxSampleCount, handler.MaxMapQ, handler.MaxReadQ,
                                    handler.AverageDirt, handler.AverageMapQ, handler.AverageReadQ,
                                    handler.AverageDirtGood, handler.AverageMapQGood, handler.AverageReadQGood);
                            Console.WriteLine("---------------------"+handler.AverageDirtGood + "--------------------");
                        }
                    }
                    else
                    {
                        // do not parse any more sequences, do not collect $500
                        Console.WriteLine("sequence loop aborted, will not parse any more");
                        break;
                    }
            }
            Dispatcher.BeginInvoke(
                System.Windows.Threading.DispatcherPriority.Normal,
                new QuickDelegate(UpdateGui_FinishedParsing));

            // Close handler and return
            Console.WriteLine("finished, no more seqs");
            handler.SetComplete(); // this should block until ok to finish
            // there might be a background thread still doing stuff
            handler.Dispose();
            parser.Dispose();
            
            return null;
        }

        private Queue<SAMAlignedSequence> inputQueue;

        private System.Delegate UpdateProgress(int index)
        {
            if(index >= ProgressBar.Minimum && index <= ProgressBar.Maximum)
            {
                ProgressBar.Value = index;
                double percent = Math.Round((index / (double)ProgressBar.Maximum * 100), 2);
                ProgressBar.ToolTip = "Analysed " + index.ToString() + " out of " + ProgressBar.Maximum.ToString() + " sequences ("+percent.ToString() + "%)";
                //Console.WriteLine("**********TT:"+ProgressBar.ToolTip);
                LoadingBarLabel.Content = ProgressBar.ToolTip;
            }
            else
            {
                throw new Exception("Invalid progress bar value: "+ProgressBar.Minimum + " : " + ProgressBar.Maximum);
            }
            if (index == ProgressBar.Maximum)
            {
                ProgressBar.Visibility = System.Windows.Visibility.Hidden;
                AbortButton.Visibility = System.Windows.Visibility.Hidden;
                LoadingBarLabel.Visibility = System.Windows.Visibility.Hidden;
            }
            return null; 
        }

        public System.Delegate DrawChart(double[] data)
        {
            SolidColorBrush myTransparentBrush = new SolidColorBrush();
            myTransparentBrush.Color = (Color)ColorConverter.ConvertFromString("Transparent");

            // Chart frame
            Chart c = new Chart();
            c.Name = "PieChart";
            c.VerticalAlignment = VerticalAlignment.Top;
            c.Margin = new System.Windows.Thickness(0, 0, 0, 0);
            c.Height = 200;
            c.Width = 200;
            c.Foreground = myTransparentBrush;
            c.BorderBrush = myTransparentBrush;

            // Chart data
            PieSeries s = new PieSeries();
            s.Title = "Frequencies";
            s.Margin = new System.Windows.Thickness(0, -40, 0, -20);
            c.Foreground = myTransparentBrush;
            c.Background = myTransparentBrush;
            c.BorderBrush = myTransparentBrush;
            s.IndependentValueBinding = new Binding("Key");
            s.DependentValueBinding = new Binding("Value");
            c.Series.Add(s);

            // Legend style
            Style style = new Style(typeof(Control));
            style.Setters.Add(new Setter(Chart.HeightProperty, new Binding("0")));
            style.Setters.Add(new Setter(Chart.WidthProperty, new Binding("0")));
            style.Setters.Add(new Setter(Chart.BorderBrushProperty, new Binding("Transparent")));
            style.Setters.Add(new Setter(Chart.BackgroundProperty, new Binding("Transparent")));
            c.LegendStyle = style;

            // Bind data
            KeyValuePair<string, double>[] dataSeries = new KeyValuePair<string, double>[data.Length];
            int i = 0;
            foreach (double d in data)
            {
                dataSeries[i] = new KeyValuePair<string, double>("Ploidy " + (++i), d);
            }
            ((PieSeries)c.Series[0]).ItemsSource = dataSeries;

            // Add chart
            wpMain.Children.Add(c);
            return null;
        }


        private void Window_Closing(object sender, System.ComponentModel.CancelEventArgs e)
        {
            if(isProcessingFile)
            {
                MessageBoxResult result = MessageBox.Show("Data is currently being output to file. \nIf you close this window the output file will be incomplete.\n\nAre you sure you want to close this window?",
                "Warning", MessageBoxButton.YesNo, MessageBoxImage.Warning);
                if (result == MessageBoxResult.No)
                {
                    e.Cancel = true;
                }
            }
        }

        private void NumericOnly_TextChanged(object sender, TextChangedEventArgs e)
        {
            // todo fixme - multiple ways of doing this, this is probably not the most efficient
            TextBox t = sender as TextBox;
            string s = "";
            int dec = 0;
            foreach (Char c in t.Text.ToCharArray())
            {
                s += Char.IsDigit(c) || (c == '.' && dec == 0)? c.ToString() : "";
                if (c == '.')
                {
                    dec++;
                }
            }
            t.Text = s;
            
        }       

    }
}
