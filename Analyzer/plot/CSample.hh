#ifndef CSAMPLE_HH
#define CSAMPLE_HH

#include <string>

//
// helper class to handle sample inputs
//
class CSample
{
public:
  CSample(){}
  CSample(const std::string mylabel, const int fillcol, const int linecol):label(mylabel),fillcolor(fillcol),linecolor(linecol){}
  ~CSample(){}
  
  std::string label;           // plot item label
  int fillcolor;               // plot item color
  int linecolor;               // plot item line color
  vector<std::string> fnamev;  // ntuple files
};
#endif
