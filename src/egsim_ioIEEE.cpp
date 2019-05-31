using namespace std;

#include "egsim.hpp"

#define MAX_CHARS_PER_LINE 200
//#define IEEE_STRICT
//#define VERBOSE
//#define DEBUG

/* 

  I/O function based on a file formats with
  - all bus variables are given at once
  - all line variables are given at once

*/
void getTextFile(const char *filein, vector<string> &ftext) ;
void getSection(vector<string> model, string sname, int *sst, int *sen) ;

void egsim::datainputIEEE(const char *filename, Array1D<bus> &bdata, Array1D<line> &ldata)
{
  vector<string> modeltext;
  getTextFile(filename, modeltext);

#ifdef DEBUG
  cout<<"--------------Input File------------------------------"<<endl;
  for ( int i=0; i<modeltext.size(); i++)
    cout<<i<<":"<<modeltext[i]<<endl;
#endif

  /* Determine start/end of bus and line data */
  int BusStart, BusEnd, LineStart, LineEnd;
  getSection(modeltext,string("BUS"),   &BusStart, &BusEnd );
  getSection(modeltext,string("BRANCH"),&LineStart,&LineEnd);
 
#ifdef VERBOSE
  cout<<"--------------Bus/Line starts/ends--------------------"<<endl;
  cout<<"Bus data :  "<<BusStart<<"-"<<BusEnd<<endl;
  cout<<"Line data: "<<LineStart<<"-"<<LineEnd<<endl;
  cout<<"------------------------------------------------------"<<endl;
#endif    

  /* power rating */
  double prat=atof(modeltext[0].substr(31,6).c_str()) ;

  /* Parsing bus data */
  int nbus=BusEnd-BusStart+1;
#ifdef VERBOSE
  printf( "-----------------------------------------------\n");
  cout << " Power rating: "<<prat<<endl;
  printf( " Parsing bus data for %d buses\n",nbus);
#endif    

  bdata.Resize(nbus);
  for ( int i = 0; i<nbus; i++ ) {
    /* Read data for bus "i" */
#ifdef IEEE_STRICT
    int    busno   = atoi(modeltext[BusStart+i].substr( 0, 4).c_str()) ; /* Bus number        */
    string busname = modeltext[BusStart+i].substr(6,11) ;                /* Bus name          */
    int    ano     = atoi(modeltext[BusStart+i].substr(18, 2).c_str()) ; /* Area number       */
    int    lzno    = atoi(modeltext[BusStart+i].substr(20, 3).c_str()) ; /* Loss zone number  */
    int    btype   = atoi(modeltext[BusStart+i].substr(24, 2).c_str()) ; /* Bus type          */
    double vmag    = atof(modeltext[BusStart+i].substr(27, 6).c_str()) ; /* Voltage magnitude */
    double vang    = atof(modeltext[BusStart+i].substr(33, 7).c_str()) ; /* Voltage angle     */
    double pload   = atof(modeltext[BusStart+i].substr(40, 9).c_str()) ; /* p load            */
    double qload   = atof(modeltext[BusStart+i].substr(49,10).c_str()) ; /* q load            */
    double pgen    = atof(modeltext[BusStart+i].substr(59, 8).c_str()) ; /* p gen             */
    double qgen    = atof(modeltext[BusStart+i].substr(67, 8).c_str()) ; /* q gen             */
    double qmax    = atof(modeltext[BusStart+i].substr(90, 8).c_str()) ; /* q max             */
    double qmin    = atof(modeltext[BusStart+i].substr(98, 8).c_str()) ; /* q min             */
    double Gshunt  = atof(modeltext[BusStart+i].substr(106,8).c_str()) ; /* Gshunt            */
    double Bshunt  = atof(modeltext[BusStart+i].substr(114,8).c_str()) ; /* Bshunt            */
#else
    std::stringstream parseDbl(modeltext[BusStart+i]);
    int    busno  ; parseDbl >> busno  ; /* Bus number        */
    string busname; parseDbl >> busname; /* Bus name          */
    string bustmp ; parseDbl >> bustmp ; /* Dummy var         */
    int    ano    ; parseDbl >> ano    ; /* Area number       */
    int    lzno   ; parseDbl >> lzno   ; /* Zone number       */
    int    btype  ; parseDbl >> btype  ; /* Bus type          */
    double vmag   ; parseDbl >> vmag   ; /* Voltage magnitude */
    double vang   ; parseDbl >> vang   ; /* Voltage angle     */
    double pload  ; parseDbl >> pload  ; /* p load            */
    double qload  ; parseDbl >> qload  ; /* q load            */
    double pgen   ; parseDbl >> pgen   ; /* p gen             */
    double qgen   ; parseDbl >> qgen   ; /* q gen             */
    double qmax   ; parseDbl >> qmax   ; /* q max             */
    double qmin   ; parseDbl >> qmin   ; /* q min             */
    double Gshunt ; parseDbl >> Gshunt ; /* Gshunt            */
    double Bshunt ; parseDbl >> Bshunt ; /* Bshunt            */
#endif

    busno--;

    /* normalize */
    vang   = vang * PI / 180.0; 
    pload /= prat; qload /= prat;
    pgen  /= prat; qgen  /= prat;
    qmax  /= prat; qmin  /= prat;

    if (btype == 3) 
      btype = SW;
    else if ( btype == 2 )
      btype = PVgn;
    else if ( btype == 0 )
      btype = PQld;
    else {
      cout<<"egsim::datainputIEEE() ERROR: Bus type "<<btype<<" not yet implemented"<<endl;
      exit(1);
    }

    pstruct vma; vma.m = vmag; vma.a = vang;
    complex<double> pqgen(pgen,qgen), pqload(pload,qload), GBs(Gshunt,Bshunt) ; 
    bdata(busno).setup(busno, vma, pqgen, pqload, GBs, btype, qmax, qmin);

#ifdef VERBOSE
    printf("%4d %.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %2d %7.4f %7.4f\n",busno,vma.m,vma.a,pqgen,pqload,GBs,btype,qmax,qmin) ;
#endif    

  } /* done reading bus data */

  /* Reading line data */
#ifdef VERBOSE
  printf( "-----------------------------------------------\n");
#endif    
  int nlin=LineEnd-LineStart+1;;
  ldata.Resize(nlin);
  for ( int i = 0; i<nlin; i++ ) {
#ifdef IEEE_STRICT
    /* Read data for line "i" */
    int    frombus = atoi(modeltext[LineStart+i].substr( 0, 4).c_str()) ; /* From Bus          */
    int    tobus   = atoi(modeltext[LineStart+i].substr( 5, 4).c_str()) ; /* To Bus            */
    int    ano     = atoi(modeltext[LineStart+i].substr(10, 2).c_str()) ; /* Area number       */
    int    lzno    = atoi(modeltext[LineStart+i].substr(12, 3).c_str()) ; /* Loss zone number  */
    int    lcirc   = atoi(modeltext[LineStart+i].substr(16, 1).c_str()) ; /* Circuit           */
    int    ltype   = atoi(modeltext[LineStart+i].substr(18, 1).c_str()) ; /* Line type         */
    double r       = atof(modeltext[LineStart+i].substr(19,10).c_str()) ; /* Resistance        */
    double x       = atof(modeltext[LineStart+i].substr(29,11).c_str()) ; /* Reactance         */
    double chrg    = atof(modeltext[LineStart+i].substr(40,10).c_str()) ; /* Line charging     */
    double tap     = atof(modeltext[LineStart+i].substr(76, 6).c_str()) ; /* Tap ratio         */
    double psa     = atof(modeltext[LineStart+i].substr(83, 7).c_str()) ; /* Phase shift angle */
#else
    std::stringstream parseDbl(modeltext[LineStart+i]);
    /* Read data for line "i" */
    int    frombus;  parseDbl>>frombus  ; /* From Bus          */
    int    tobus  ;  parseDbl>>tobus    ; /* To Bus            */
    int    ano    ;  parseDbl>>ano      ; /* Area number       */
    int    lzno   ;  parseDbl>>lzno  ; /* Loss zone number  */
    int    lcirc  ;  parseDbl>>lcirc  ; /* Circuit           */
    int    ltype  ;  parseDbl>>ltype ; /* Line type         */
    double r    ; parseDbl>>r    ; /* Resistance        */
    double x    ; parseDbl>>x    ; /* Reactance         */
    double chrg ; parseDbl>>chrg ; /* Line charging     */
    int    ltmp1; parseDbl>>ltmp1;
    int    ltmp2; parseDbl>>ltmp2;
    int    ltmp3; parseDbl>>ltmp3;
    int    ltmp4; parseDbl>>ltmp4;
    int    ltmp5; parseDbl>>ltmp5;
    double tap  ; parseDbl>>tap  ; /* Tap ratio         */
    double psa  ; parseDbl>>psa  ; /* Phase shift angle */
#endif

    psa *= PI / 180.0;
    frombus--; tobus--;
    ldata(i).setup(frombus,tobus,complex<double>(r,x),chrg,tap,psa);
    ldata(i).set_on() ;

#ifdef VERBOSE
    printf( "%4d %4d %f %f %f %f %f\n",frombus,tobus,r,x,chrg,tap,psa);
#endif

  } /* done reading line data */

  return ;

}

void getTextFile(const char *filein, vector<string> &ftext) {

  ifstream fin ;
  fin.open(filein); 
  if (!fin.good()) {
    cout << "getTextFile() Could not open "<<filein<<endl<<flush ;
    exit(0);
  }

  while (!fin.eof()) {
    // read an entire line into memory
    char buf[MAX_CHARS_PER_LINE];
    fin.getline(buf, MAX_CHARS_PER_LINE);
    string bufstr(buf) ;
    ftext.push_back(bufstr);
  }
  fin.close();
  return ;

}

void getSection(vector<string> model, string sname, int *sst, int *sen) {
  
  bool foundSection = false;
  int i = 0 ;
  *sst = *sen = 0;
  while ( ( i < (int) model.size() ) && ( *sen == 0 ) ){
    size_t sFound=model[i].find(sname) ;
    if ( sFound != string::npos) {
      foundSection = true ;
      *sst = i+1;
    }
    else if ( foundSection ) {
      size_t eFound=model[i].find(string("-999")) ;
      if ( eFound != string::npos) {
        foundSection = false ;
        *sen = i-1;
      }
    }
    i++ ;
  }

  return ;
  
}
