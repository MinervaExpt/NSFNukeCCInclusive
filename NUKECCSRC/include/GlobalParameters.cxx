#ifndef GLOBALPARAMETERS_cxx
#define GLOBALPARAMETERS_cxx 1

#include "include/GlobalParameters.h"
#include "include/GlobalIncludes.h" 

using namespace CCQENU_ANA; 

//--Default Constructor 
GlobalParameters::GlobalParameters() 
{

  m_useFluxConstraint = false; 
 //  m_neutrinoMode = true; 
  m_usePPFX1Flux  = true;  
  m_analysisType = k2DAntiNuLepton;
// m_analysisType = k2DNuLepton;
/*
  std::string mnvRelease = getenv("MINERVA_RELEASE");
  //std::cout << "Release no " << mnvRelease << std::endl; 

  std::size_t found = mnvRelease.find("r"); 
  found = mnvRelease.find_first_of( "0123456789", found+1); 
  
  stringstream ss; 
  string str2;
  ss << mnvRelease[found]; 
  ss >> str2; 
  
  // RESURRECTION processing was frozen at a v10r7pX release, EROICA started on a v10r8pX 
  if( std::atoi(str2.c_str()) < 8 ) { 
    m_single_nu_mc = true; 
  } else { 
    m_single_nu_mc = false; 
  }
*/
}

GlobalParameters& GlobalParameters::Get() 
{ 
  
  static GlobalParameters singleton; 
  return singleton; 
}

void GlobalParameters::Print()
{ 

  std::cout << std::boolalpha << m_useFluxConstraint << std::endl; 
  //  std::cout << std::boolalpha << m_neutrinoMode << std::endl; 
  
}

GlobalParameters& globalParameters() 
{

  return GlobalParameters::Get(); 
}


#endif 
