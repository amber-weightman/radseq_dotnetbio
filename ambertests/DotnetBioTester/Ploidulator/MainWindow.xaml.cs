using Bio.Algorithms.Metric;
using Bio.IO.BAM;
using Bio.IO.SAM;
using System;
using System.Collections.Generic;
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

        public delegate System.Delegate MyDelegate(string a, string b, string c, string d, string e, string f, string g);
        public delegate System.Delegate OneArgDelegate(string message);
        public delegate System.Delegate SimpleDelegate(ref ClusterMetricHandlerPloidulator handler);

        public MainWindow()
        {
            InitializeComponent();
        }

     

        // On click of GO button, begin reading file
        private void Button_Click_Go(object sender, RoutedEventArgs e)
        {
            if (DataFileATextbox.Text != null)
            {
                MyDelegate handler = ParseBAMMetric;
                handler.BeginInvoke(DataFileATextbox.Text, ExpectedPloidyTextbox.Text, 
                    DirtCutoffTextbox.Text, AlignmentQualCutoffTextbox.Text,
                    ReadQualCutoffTextbox.Text, NumSamplesTextbox.Text, NumClustersToParseTextbox.Text, null, null);
            }
        }


        // Demo that window thread is live
        private void Button_Click_Circle(object sender, RoutedEventArgs e)
        {
            Dispatcher.BeginInvoke(
                System.Windows.Threading.DispatcherPriority.Normal,
                new OneArgDelegate(DrawCircle),
                "");
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
            } 
        }

        private System.Delegate UpdateGui(string message)
        {
            DataFileATextbox.Text = "";
            FileExpander.IsExpanded = false;
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
        private System.Delegate LoadMorePies(ref ClusterMetricHandlerPloidulator handler)
        {
            if(handler != null)
            {
                handler.MaxNumCharts = handler.MaxNumCharts + 10;
            }
            
            return null;
        }

        private System.Delegate ParseBAMMetric(string filename, string ploidy, string dirtCutoff, string alignQualCutoff, 
            string readQualCutoff, string numSamples, string numClustersToParse)
        {
            string newName = filename.Split(new char[]{'.'})[0];
            handler = new ClusterMetricHandlerPloidulator(newName + "_pl", Convert.ToInt32(ploidy), Convert.ToDouble(dirtCutoff), 
                Convert.ToInt32(alignQualCutoff), readQualCutoff, Convert.ToInt32(numSamples), Convert.ToInt32(numClustersToParse), Dispatcher, wpMain); 
            handler.WriteToFile(false); // todo oops, it still opens file
            
            BAMParser parser = new BAMParser(handler, false);

            Dispatcher.BeginInvoke(
                System.Windows.Threading.DispatcherPriority.Normal,
                new OneArgDelegate(UpdateGui),
                "");

            parser.Parse(filename);
            return null;
        }

        ClusterMetricHandlerPloidulator handler = null;

        private void Button_Click_LoadMore(object sender, RoutedEventArgs e)
        {
            if (this.handler != null)
            {
                SimpleDelegate del = LoadMorePies;
                del.BeginInvoke(ref handler, null, null);
            }
            
        }


    }
}
