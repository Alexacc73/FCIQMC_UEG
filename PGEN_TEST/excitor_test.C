#include "planeWaves.H"     /* Includes <iostream>  <set>  <math.h>  <stdlib.h>  <algorithm> */
#include "UEGHamiltonian.H" /* includes <bitset>  <climits>  <boost/math/special_functions/binomial.hpp>  <vector> */

#include <fstream>
#include <time.h>
#include <set>
#include <map>


/*
*-----> Parameters for initial setup of UEG <----- 
*/

/** A number, nobody knows what it means */
const double PI = 3.141592653589793;

/** rs controls density of the Electron Gas. 
* "rs" is the radius of the sphere whose volume is that of the cell divided by the number of electrons*/
const double rs = 0.5; 

/** Total number of electrons -> half allocated Alpha spin, the other half allocated Beta Spin.*/
const double numElectrons = 4; 
const int INTelectrons = numElectrons ;

/** Kc_CUTTOFF is the kinetic energy cutoff for the plane wave basis orbitals.
E.g, a cutoff of "2" will allow the orbital [4 0 0] but not [5 0 0]. Set cutoff = 2.4 for 57 Orbitals (114 Spin Orbitals) */
const double Kc_CUTTOFF = 2 ; 



/*
 *-----> Parameters for Main FCIQMC Algorithm <----- 
 */

/** delt is the Imaginary timestep for the propogation of the "walker" population */
const double delt = 0.00025 ;

/** Zeta is a damping parameter which controls the agressiveness of the "shift" in the variable shift mode of the algorithm */
const double zeta = 0.005 ;

/** AShift controls how frequently the shift is changed in response to the population in the variable shift mode (AShift = 1 means every step) */
const int AShift = 5 ;

/** Number of steps after which to terminate the algorithm*/
const int numSteps = 1000000;

/** After "walker critical" walkers have been spawned after a complete cycle (post annihilation) the variable shift mode is turned on */
const int walkerCritical = 100000;

long int pow2Array [ORB_SIZE];


//--------------------------------------------------
const bool ParaSPINexcite = false;
//--------------------------------------------------



/**
Simply raises the input argument to the power 2, thus returning \f$ 2^x \f$
*/
inline constexpr std::uint64_t INLpow2 (std::uint64_t x)
{
    return std::uint64_t(1) << x;
}


/**
* This function purely resturns the position of a unique determinant within the list of Alpha and Beta determinants.
* It is necessary to know this, in order that we can add a new walker to the correct determinant (Our unique list is not
* sorted according to the walker number list, but the Alpha and Beta lists ARE sorted to associate with the walker number list)
*/
inline  int INLgetPositionInList(std::pair<long int, long int> uniqueDet, std::vector<long int> alphaDetList, std::vector<long int> betaDetList)
{
    int DETLIST_size = alphaDetList.size();
    int pos ;
    int j = 0;
    bool not_found = true ;

    long int alphaBIN = uniqueDet.first ;
    long int betaBIN = uniqueDet.second ;
    while( (j<DETLIST_size) && (not_found == true)  ){
        if( (alphaBIN == alphaDetList[j]) && (betaBIN == betaDetList[j]) ){
            pos = j ;
            not_found = false ;
        }
        j++ ;
    }
    if(not_found == true){
        std::cout << "OOPS POTISION FIND WENT WRONG!!!!!" << std::endl;
    }
    return pos ;
}


/**
* This function simply returns the number of bits, i.e, the number of ones in a binary number.
*/
inline size_t oneBitCount(long int n){
    std::bitset<sizeof(size_t) * CHAR_BIT> b(n);
    return b.count();
}


/**
* This inline function simply converts a base 10 integer (which represents the unique determinant)
* into its binary representation, but with type string. 
*/
inline  std::string INLdecimalToBinary(long int decimal)
{
    std::string binaryNum;
    binaryNum = std::bitset<ORB_SIZE>(decimal).to_string() ;
    return binaryNum;
}

/**
* This inline function simply converts a string (of length defined at compile time) which represents the unique determinant
* into its decimal representation, type long int (up to 2^63).
*/
inline long int INLbinaryToDecimal(std::string binaryString)
{
    long int decimal;
    decimal = std::bitset<ORB_SIZE>(binaryString).to_ulong() ;
    return decimal;
}


int totalWalkerNumber( std::vector<int>& trueWalkerList){
  int count = 0;
  int numDets = trueWalkerList.size() ;
  for(int det = 0; det < numDets; det++){
    count += abs(trueWalkerList[det]);
  }
  return count;
}










int main(void){
    srand(492831);

    for(int i = 0; i<ORB_SIZE; i++){
        pow2Array[i] = INLpow2(i);
    }

	const double cellVolume = (numElectrons*PI*4.0*rs*rs*rs)/3.0;
	const double cellLength = getCellLength(cellVolume); /* Return cell length in units of rs */
	int spinOrbitals;
    double kpoints[500][3];

    createPlaneWaveOrbitals(kpoints, Kc_CUTTOFF, spinOrbitals);

    std::cout << "Cell length = " << cellLength << '\n' << std::endl;
    std::cout << "Number spin orbitals = " << spinOrbitals << std::endl;
    std::cout << "Number of orbitals = " << ORB_SIZE  << std::endl;
    std::cout << "Number of Electrons = " << numElectrons << '\n' << std::endl;
    if(spinOrbitals/2 != ORB_SIZE){
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::cout << "ERROR: Orbital size is not correct: Alter in file < GLOBAL_ORBITAL.H > " << std::endl;
        std::cout << "---> SET < ORB_SIZE > TO BE : " << spinOrbitals/2 << std::endl;
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        return 1;
    }

    double KEsortedKpoints[ORB_SIZE][3];
    /* Create a set of UNIQUE orbitals, in ascending Kinetic energy order */
    kPointsEnergySort(kpoints, KEsortedKpoints, ORB_SIZE);
    std::cout << "ORDERED LIST: " << std::endl;
    PRINTORBITALS(ORB_SIZE, KEsortedKpoints);

    /* -----> DEFINE IMPORTANT LISTS <----- */
    std::set< std::pair<long int, long int> > uniqueDeterminantSet ;
    std::pair< std::set< std::pair<long int, long int> >::iterator , bool > result;

    std::vector<long int> alphaDetsBinary;
    std::vector<long int> betaDetsBinary;

    std::vector<int> walkerList;
  

    long int zero = 0;
    std::string HF_STRING = INLdecimalToBinary(zero) ; 
    for(int el = 0; el < INTelectrons/2; el++){
        HF_STRING[ORB_SIZE-1-el] = '1';
    }

    long int HF_binary = INLbinaryToDecimal(HF_STRING);
    result = uniqueDeterminantSet.insert( std::make_pair(HF_binary, HF_binary) ) ; 
    if(result.second == true){
        std::cout << "INITIAL INPUT SUCCESSFUL" << std::endl;
        std::pair<long int, long int> iter = *result.first;
        std::cout << "HF Alpha + Beta DETERMINANT = " << iter.first << " " << iter.second << std::endl;
        //alphaDetsBinary.push_back( HF_binary ) ;
        //betaDetsBinary.push_back( HF_binary ) ;
    }
    double EnergyHFnonConst;
    Di_H_Di(cellLength, INTelectrons, HF_binary, HF_binary, KEsortedKpoints, EnergyHFnonConst ) ;
    const double HFEnergy = EnergyHFnonConst;
    std::cout << " < D_I | H | D_I > Element for HF = " << HFEnergy << '\n' << std::endl;
    uniqueDeterminantSet.clear();




    int TESTS = 1000000;    
    for(int i = 0; i<TESTS; i++){
        int index_i = 0, index_a = 0, index_j = 0, index_b = 0;
        int sign = 0;
        int position;
        std::pair<long int, long int> iter;
        long int alphaDet = HF_binary;
        long int betaDet = HF_binary;

        if(ParaSPINexcite){
            excitationSameSpinij_ab(index_i, index_a, index_j, index_b, INTelectrons, ORB_SIZE, alphaDet, KEsortedKpoints, sign) ;
        }
        else{
            excitationAlpha_iaBeta_jb(index_i, index_a, index_j, index_b, INTelectrons, ORB_SIZE, alphaDet, betaDet, KEsortedKpoints, sign) ;
        }

        //excitationAlpha_iaBeta_jb(index_i, index_a, index_j, index_b, INTelectrons, ORB_SIZE, HF_binary, HF_binary, KEsortedKpoints, sign) ;
        if(index_b != -1){
            //std::cout << "s_i --> s_a = " << index_i << " --> " << index_a << std::endl;
            //std::cout << "i --> a = [" << KEsortedKpoints[index_i][0] << KEsortedKpoints[index_i][1] << KEsortedKpoints[index_i][2] ;
            //std::cout << "] --> [" << KEsortedKpoints[index_a][0] << KEsortedKpoints[index_a][1] << KEsortedKpoints[index_a][2] <<"]"<< std::endl;
            //std::cout << "s_j --> s_b = " << index_j << " --> " << index_b << std::endl;
            //std::cout << "j --> b = [" << KEsortedKpoints[index_j][0] << KEsortedKpoints[index_j][1] << KEsortedKpoints[index_j][2] ;
            //std::cout << "] --> [" << KEsortedKpoints[index_b][0] << KEsortedKpoints[index_b][1] << KEsortedKpoints[index_b][2] <<"]"<< std::endl;
            if(ParaSPINexcite){
                alphaDet -= pow2Array[ index_i] ; // = '0' ;
                alphaDet += pow2Array[ index_a] ; // = '1' ;
                alphaDet -= pow2Array[ index_j] ; // = '0' ;
                alphaDet += pow2Array[ index_b] ; // = '1' ;
            }
            else{
                alphaDet -= pow2Array[ index_i] ; // = '0' ;
                alphaDet += pow2Array[ index_a] ; // = '1' ;
                betaDet -= pow2Array[ index_j] ; // = '0' ;
                betaDet += pow2Array[ index_b] ; // = '1' ;

            }
            //std::cout << "new Alpha Det = " << alphaDet << "| in string form = " << INLdecimalToBinary(alphaDet) <<  std::endl;
            //std::cout << "new Beta Det = " <<  betaDet << "| in string form = " << INLdecimalToBinary(betaDet) << '\n' << std::endl;

            result = uniqueDeterminantSet.insert( std::make_pair(alphaDet, betaDet) ) ;
            if(result.second == true){/* New Determinant Found */
                walkerList.push_back(1) ;
                alphaDetsBinary.push_back(alphaDet) ;
                betaDetsBinary .push_back(betaDet) ;
            }
            if(result.second == false){ /* existing determinant found */
                iter = *result.first ;
                /*Find position of determinant in list*/
                position = INLgetPositionInList(iter, alphaDetsBinary, betaDetsBinary) ;
                walkerList[position] += 1;

            }
        }
    }
    double pGen;
    //pGen = (1.0/((numElectrons/2.0)*(numElectrons/2.0)) ) * ( 1.0/( float(ORB_SIZE) - (numElectrons/2)) ) ;

    if(ParaSPINexcite){
        pGen = (2.0/( (numElectrons/2.0) *((numElectrons/2.0)-1.0)) ) * ( 2.0/(float(ORB_SIZE) - (numElectrons/2)) ) ;
    }
    else{
        pGen = (1.0/((numElectrons/2.0)*(numElectrons/2.0)) ) * ( 1.0/( float(ORB_SIZE) - (numElectrons/2)) ) ;
    }

    //pGen = pGen/2.0 ;
    

    std::cout << "Number of double excitations from HF det  = " << uniqueDeterminantSet.size() << std::endl;
    std::cout << "walker list size  = " << walkerList.size() << std::endl;
    std::cout << "Alpha det size  = " << alphaDetsBinary.size() << std::endl;
    std::cout << "Total number of valid doubles produced = " << totalWalkerNumber(walkerList) << std::endl;

    double sum = 0;
    float counter = 0;
    for (std::vector<int>::const_iterator i = walkerList.begin(); i != walkerList.end(); ++i){
        std::cout << "Number generated = "<<  *i << '\n';
        sum += *i ;
        counter += 1 ;
    }
    std::cout << " " << std::endl;
    std::cout << "first Alpha det = " << INLdecimalToBinary(alphaDetsBinary[0]) <<  std::endl;
    std::cout << "first beta det = " << INLdecimalToBinary(betaDetsBinary[0]) <<  std::endl;

    std::cout << '\n' << "Number det 1 created / total number of attempts = " << float(walkerList[0]) / float(TESTS) << std::endl;
    std::cout << "--------------------------" << std::endl;
    std::cout << "Deterministic pGen = " << pGen << std::endl;
    std::cout << "Empirical pGen = " << sum/counter / float(TESTS) << std::endl;


    return 1;

}