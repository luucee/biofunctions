viterbi <- function(transProbs,emissionProbs,States,startProbs, observation)
{
    transProbs[is.na(transProbs)] = 0
    emissionProbs[is.na(emissionProbs)] = 0
    nObservations = length(observation)
    nStates = length(States)
    v = array(NA, c(nStates, nObservations))
    p = array(NA, c(nStates, nObservations))
    dimnames(v) = list(states = States, index = 1:nObservations)
    dimnames(p) = list(states = States, index = 1:nObservations)
    v[States, 1] = log(startProbs[States] * emissionProbs[States,as.character(observation[1])])
    p[States, 1] = 0
    for (k in 2:nObservations) {
        for (state in States) {
            #cat(k,state,"\n")
            temp = v[States, k - 1]+log(transProbs[States,state])
            maxi=which.max(temp)
            v[state, k] = log(emissionProbs[state, as.character(observation[k])]) + temp[maxi]
            p[state, k] = maxi
        }
        if (sum(v[States, k]==-Inf) == nStates) {
            cat("impossible to find an optimal transition for the observation: ",observation[k-1],"-->",observation[k],"\n",file=stderr())
            #caso di transizione non riconosciuta nel training set
            #ricomputo senza considerare la transizione... sperando in una buona sorte!
            temp = v[States, k - 1]
            maxi=which.max(temp)
            for (state in States) {
                v[state, k] = log(emissionProbs[state, as.character(observation[k])]) + temp[maxi]
                p[state, k] = maxi
            }
        }
    }
    viterbiPath = rep(NA, nObservations)
    state.k<-which.max(v[States, nObservations])
    viterbiPath[nObservations]<-States[state.k]
    for (k in nObservations:2) {
        state.k<-p[state.k,k]
        viterbiPath[k-1] =States[state.k]
    }
    return(viterbiPath)
}

train <- function(ob,st) {
  S<-unique(st)
  SY<-unique(ob)
  M<-matrix(0,nrow=length(S),ncol=length(S))
  E<-matrix(0,nrow=length(S),ncol=length(SY))
  rownames(M)<-S
  colnames(M)<-S
  rownames(E)<-S
  colnames(E)<-SY
  
  for(i in 2:length(ob)) {
    M[st[i-1],st[i]]<-M[st[i-1],st[i]]+1
    E[st[i],ob[i]]<-E[st[i],ob[i]]+1
  }
  sumM<-apply(M,1,sum)
  M<-M/sumM
  sumE<-apply(E,1,sum)
  E<-E/sumE
  return(list(M=M,E=E))
}
