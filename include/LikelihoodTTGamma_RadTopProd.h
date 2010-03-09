/*!
 * \class KLFitter::LikelihoodTTGamma_RadTopProd
 * \brief A class implementing a likelihood for radiative top production
 * \author Johannes Erdmann
 * \version 1.3
 * \date 26.02.2010
 *
 * This class represents a likelihood for radiative top production
 */

// --------------------------------------------------------- 

#ifndef __LIKELIHOODTTGAMMA_RADTOPPROD__H
#define __LIKELIHOODTTGAMMA_RADTOPPROD__H

// --------------------------------------------------------- 

#include "LikelihoodTTGamma.h" 

// --------------------------------------------------------- 

/*!
 * \namespace KLFitter
 * \brief The KLFitter namespace
 */
namespace KLFitter
{

	class LikelihoodTTGamma_RadTopProd : public KLFitter::LikelihoodTTGamma
	{
		
	public: 
		
		/** 
		 * The default constructor. 
		 */ 
		LikelihoodTTGamma_RadTopProd(); 
		
		/**
		 * The default destructor.
		 */
		virtual ~LikelihoodTTGamma_RadTopProd(); 

	}; 

} // namespace KLFitter 

// --------------------------------------------------------- 

#endif // __LIKELIHOODTTGAMMA_RADTOPPROD__H
