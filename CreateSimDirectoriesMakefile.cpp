// CreateSimDirectoriesMakefile.cpp
#include <fstream>
using namespace std;

string ReadEntryStringFromConfigFile(int entryId) {
  fstream infile;
  infile.open("configfile",ios::in);
  string temp1,temp2;
  for (int i=0; i<entryId; i++) 
    infile>>temp1>>temp2; //move to correct entry
  infile >>temp1>>temp2;
  return temp2; 
}

void CreateSimDirectoriesMakefile() {
  string SimName=ReadEntryStringFromConfigFile(11);
  SimName="SimulationData/"+SimName;
  fstream outfile;
  outfile.open("CreateSimDirectories.make",ios::out);
  outfile<<"CreateSimDirectories:\n";
  outfile<<"\tmkdir "<<SimName<<endl;
  outfile<<"\tmkdir "<<SimName<<"/DetectData\n";
  outfile<<"\tmkdir "<<SimName<<"/MatrixData\n";
  outfile.close();
}

int main() {
  CreateSimDirectoriesMakefile();
}


