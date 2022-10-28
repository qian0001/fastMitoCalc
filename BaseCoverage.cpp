/*baseCoverage.cpp by TJ Butler
 * 8/9/13: removed X and Y from average coverage calculation
 * 11/13/13: added necessary working directory option
 * 6/17/14: Annotated code for distribution
 * 7/2/14: Changed name to baseCoverage.cpp from baseCoveragePiped.cpp,
 * did final look through of code annotation
 * 5/10/2016 Modified the code for use with single chromosome. Also changed the print out of text file by YQ
 * last modified 05/10/2016 
*/

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <vector>
#include <map>

using namespace std;

int find_autos(string chromNumber, string prefix) {
  string autos[]={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"};
  for(int i=0;i<22;i++) {
     string chrom = prefix+autos[i];
     if(chromNumber.compare(chrom)==0) {
        return 0;
     }
  }
  return 1;
}


void basecoverage(ifstream& infile,double& coverageAutosomal, double& coverageMT, int AutosomeUsed[],
       int& AutosomeCnt, int& bpcnt, string prefix, string mtGenome, int& mtbases)
//Function to determine the average base coverage of autosomal and mitochrondria chromosomes
//Reads in lines from SAMtools depth, and calculates average coverage of that chromosome and
//returns the info for that.
//This takes number of bases read into account, as N-region likely does not
//have any coverage of bases thus throw out those bases altogether.
{
	string chromNumber;
	double autobases = 0;
	double autosum = 0;
	double autoavg = 0;
	double mtsum = 0;
	double mtavg = 0;
	int    curCoverage = 0;
	double curBaseNo = 0;
        int chr;

	//Loop through SAMtools depth output file line by line
	infile >> chromNumber >> curBaseNo >> curCoverage;
        chr = atoi(chromNumber.c_str() + prefix.length());
        AutosomeUsed[0] = chr;
        AutosomeCnt = 1;

        while(infile.good()) {
          int findA = find_autos(chromNumber,prefix);
          if(findA!=0 && chromNumber.compare(mtGenome) != 0) {
             infile >> chromNumber >> curBaseNo >> curCoverage;
             chr = atoi(chromNumber.c_str() + prefix.length());
             continue;
          }

          if(findA==0) {
              coverageAutosomal += curCoverage;
              autobases++;

            if(chr != AutosomeUsed[AutosomeCnt-1]) {
              AutosomeUsed[AutosomeCnt] = chr;
              AutosomeCnt++;
            }
          }
          else if((chromNumber.compare(mtGenome) == 0)){
             coverageMT += curCoverage;
             mtbases++;
          }
          infile >> chromNumber >> curBaseNo >> curCoverage;
          chr = atoi(chromNumber.c_str() + prefix.length());
        }
        coverageAutosomal = (double)(coverageAutosomal/autobases);
        coverageMT = (double)(coverageMT/mtbases);
        bpcnt = autobases;

}

int main(int argc, char** argv){
  if ( (argc < 5) || (argc > 6) ) {
    cout << "\nusage: " << argv[0] << " </working/dir/> <bamFileLabel> <depthfile> <mt> <prefix>\n\n";
    cout << "*This program finds the average depth of each chromosome, 1-22 and MT.\n";
    cout << "<bamFileLabel> is a unique identifier for a file, to be appended\n";
    cout << "to uniquely label results files generated.\n";
    cout << "/working/dir/ is complete path to a working directory where result\n";
    cout << "files should be saved\n";
    cout << "If current directory is desired, type \"./\"\n\n";
    cout << "if autosomal sequences are 1,2,..22 don't use -p option. If chr1,..chr22, use -p chr\n";
    cout << "if mitochondria is named as MT don't use -m option. Otherwise use names like chrM, M, depending on how it is defined in bam file\n";
    cout << "*This program utilizes samtools' depth utility\n";
    return 1;
  }

  string workDir = argv[1];
  string bamLabel = argv[2];
  string infileName = argv[3];
  string mtGenome = argv[4];
  string prefix = "";
  if(argc==6)
     prefix = argv[5];
	
   //will come out as '<bamFileLabel>.txt'
   string outfileName = workDir + "/" + bamLabel + ".txt";
   ifstream infile;
   infile.open(infileName.c_str());
   ofstream outfile;
   outfile.open(outfileName.c_str());

   int AutosomeUsed[22];
      int AutosomeCnt;
      int bpcnt;

      //For entire SAMtools depth file, go through each chromosome and stop when at MT
      double autoZcoverage= 0;
      double MTcoverage = 0;
      int    mtbases = 0;

      basecoverage(infile,autoZcoverage, MTcoverage, AutosomeUsed,AutosomeCnt,bpcnt,prefix,mtGenome,mtbases);

      if(!mtbases) {
        cout << "No coverage in "<<mtGenome<<", make sure your mt chromosomes is named as " <<mtGenome<<endl;
        abort ();
      }

         
      infile.close();

      //Calculate copy number estimates for MT: using average calculated above
      double mtCopyNoAvg = 2*MTcoverage/(autoZcoverage);

      //Write everything to the outfile
      outfile << "mt_copy_number_avg\t" ;
      outfile << "mt_coverage\t" ;
      outfile << "autosomal_coverage\t";
      outfile << "actual basepairs covered by reads\t";
      outfile << "chrom_used_for_autosomal_coverage"<< endl;


      outfile << mtCopyNoAvg << "\t";
      outfile << MTcoverage << "\t";
      outfile << autoZcoverage << "\t";
      outfile << bpcnt << "\t";
      if(AutosomeCnt<22) {
         for (int i=0; i<AutosomeCnt-1; i++) {
             outfile << AutosomeUsed[i]<<",";
         }
         outfile << AutosomeUsed[AutosomeCnt-1] << endl;
      }
      else
          outfile << "1-22" << endl;
      outfile.close();
}

//End of program
