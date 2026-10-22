import numpy as np

def gradient_descent(data,model_func,initial_coeff):
    """
    intial_coeff of the form of np.array([[npop,nparams]])

    """
    
    stepsize=0.001
    mycounter=1
    dRSS=1
    lastRSS=10,
    npop,nparams= initial_coeff.shape
    current_coeff = initial_coeff
    
    perturbation = np.zeros((npop,nparams,npop,nparams))
    
    for i in range(npop):
        for j in range(nparams):
            perturbation[i,j,i,j] = 0.0001

    while ((mycounter < 5000) & (dRSS > 10**-6)):                                      
        
        myRSS = np.linalg.norm(data-model_func(current_coeff))**2

        dRSS= abs(lastRSS-myRSS)
        lastRSS = myRSS
        
        partialDs = np.array([[(np.linalg.norm(data-model_func(current_coeff+perturbation[i,j]))**2 - 
                     myRSS)/.0001 for j in range(nparams)] for i in range(npop)])
        

        
        current_coeff=abs(current_coeff-stepsize*partialDs)


        mycounter+=1
    return current_coeff