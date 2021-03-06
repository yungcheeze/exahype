#include "exahype/records/State.h"

#if !defined(TrackGridStatistics) && !defined(Parallel)
   exahype::records::State::PersistentRecords::PersistentRecords() {
      
   }
   
   
   exahype::records::State::PersistentRecords::PersistentRecords(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
   _maxRefinementLevelAllowed(maxRefinementLevelAllowed),
   _mergeMode(mergeMode),
   _sendMode(sendMode),
   _algorithmSection(algorithmSection),
   _hasRefined(hasRefined),
   _hasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration),
   _hasErased(hasErased),
   _hasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration),
   _hasChangedVertexOrCellState(hasChangedVertexOrCellState),
   _hasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration),
   _isTraversalInverted(isTraversalInverted) {
      
   }
   
   exahype::records::State::State() {
      
   }
   
   
   exahype::records::State::State(const PersistentRecords& persistentRecords):
   _persistentRecords(persistentRecords._maxRefinementLevelAllowed, persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._algorithmSection, persistentRecords._hasRefined, persistentRecords._hasTriggeredRefinementForNextIteration, persistentRecords._hasErased, persistentRecords._hasTriggeredEraseForNextIteration, persistentRecords._hasChangedVertexOrCellState, persistentRecords._hasModifiedGridInPreviousIteration, persistentRecords._isTraversalInverted) {
      
   }
   
   
   exahype::records::State::State(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
   _persistentRecords(maxRefinementLevelAllowed, mergeMode, sendMode, algorithmSection, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted) {
      
   }
   
   
   exahype::records::State::~State() { }
   
   std::string exahype::records::State::toString(const AlgorithmSection& param) {
      switch (param) {
         case TimeStepping: return "TimeStepping";
         case LimiterStatusSpreading: return "LimiterStatusSpreading";
         case MeshRefinement: return "MeshRefinement";
         case MeshRefinementOrLocalOrGlobalRecomputation: return "MeshRefinementOrLocalOrGlobalRecomputation";
         case LocalRecomputationAllSend: return "LocalRecomputationAllSend";
         case MeshRefinementOrGlobalRecomputation: return "MeshRefinementOrGlobalRecomputation";
         case MeshRefinementOrGlobalRecomputationAllSend: return "MeshRefinementOrGlobalRecomputationAllSend";
         case PredictionRerunAllSend: return "PredictionRerunAllSend";
      }
      return "undefined";
   }
   
   std::string exahype::records::State::getAlgorithmSectionMapping() {
      return "AlgorithmSection(TimeStepping=0,LimiterStatusSpreading=1,MeshRefinement=2,MeshRefinementOrLocalOrGlobalRecomputation=3,LocalRecomputationAllSend=4,MeshRefinementOrGlobalRecomputation=5,MeshRefinementOrGlobalRecomputationAllSend=6,PredictionRerunAllSend=7)";
   }
   std::string exahype::records::State::toString(const MergeMode& param) {
      switch (param) {
         case MergeNothing: return "MergeNothing";
         case BroadcastAndMergeTimeStepData: return "BroadcastAndMergeTimeStepData";
         case MergeFaceData: return "MergeFaceData";
         case DropFaceData: return "DropFaceData";
         case BroadcastAndMergeTimeStepDataAndMergeFaceData: return "BroadcastAndMergeTimeStepDataAndMergeFaceData";
         case BroadcastAndMergeTimeStepDataAndDropFaceData: return "BroadcastAndMergeTimeStepDataAndDropFaceData";
      }
      return "undefined";
   }
   
   std::string exahype::records::State::getMergeModeMapping() {
      return "MergeMode(MergeNothing=0,BroadcastAndMergeTimeStepData=1,MergeFaceData=2,DropFaceData=3,BroadcastAndMergeTimeStepDataAndMergeFaceData=4,BroadcastAndMergeTimeStepDataAndDropFaceData=5)";
   }
   std::string exahype::records::State::toString(const SendMode& param) {
      switch (param) {
         case SendNothing: return "SendNothing";
         case ReduceAndMergeTimeStepData: return "ReduceAndMergeTimeStepData";
         case SendFaceData: return "SendFaceData";
         case ReduceAndMergeTimeStepDataAndSendFaceData: return "ReduceAndMergeTimeStepDataAndSendFaceData";
      }
      return "undefined";
   }
   
   std::string exahype::records::State::getSendModeMapping() {
      return "SendMode(SendNothing=0,ReduceAndMergeTimeStepData=1,SendFaceData=2,ReduceAndMergeTimeStepDataAndSendFaceData=3)";
   }
   
   
   std::string exahype::records::State::toString() const {
      std::ostringstream stringstr;
      toString(stringstr);
      return stringstr.str();
   }
   
   void exahype::records::State::toString (std::ostream& out) const {
      out << "("; 
      out << "maxRefinementLevelAllowed:" << getMaxRefinementLevelAllowed();
      out << ",";
      out << "mergeMode:" << toString(getMergeMode());
      out << ",";
      out << "sendMode:" << toString(getSendMode());
      out << ",";
      out << "_algorithmSection:" << toString(getAlgorithmSection());
      out << ",";
      out << "hasRefined:" << getHasRefined();
      out << ",";
      out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
      out << ",";
      out << "hasErased:" << getHasErased();
      out << ",";
      out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
      out << ",";
      out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
      out << ",";
      out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
      out << ",";
      out << "isTraversalInverted:" << getIsTraversalInverted();
      out <<  ")";
   }
   
   
   exahype::records::State::PersistentRecords exahype::records::State::getPersistentRecords() const {
      return _persistentRecords;
   }
   
   exahype::records::StatePacked exahype::records::State::convert() const{
      return StatePacked(
         getMaxRefinementLevelAllowed(),
         getMergeMode(),
         getSendMode(),
         getAlgorithmSection(),
         getHasRefined(),
         getHasTriggeredRefinementForNextIteration(),
         getHasErased(),
         getHasTriggeredEraseForNextIteration(),
         getHasChangedVertexOrCellState(),
         getHasModifiedGridInPreviousIteration(),
         getIsTraversalInverted()
      );
   }
   
   #ifdef Parallel
      tarch::logging::Log exahype::records::State::_log( "exahype::records::State" );
      
      MPI_Datatype exahype::records::State::Datatype = 0;
      MPI_Datatype exahype::records::State::FullDatatype = 0;
      
      
      void exahype::records::State::initDatatype() {
         {
            State dummyState[2];
            
            #ifdef MPI2
            const int Attributes = 11;
            #else
            const int Attributes = 12;
            #endif
            MPI_Datatype subtypes[Attributes] = {
                 MPI_INT		 //maxRefinementLevelAllowed
               , MPI_INT		 //mergeMode
               , MPI_INT		 //sendMode
               , MPI_INT		 //_algorithmSection
               , MPI_CXX_BOOL		 //hasRefined
               , MPI_CXX_BOOL		 //hasTriggeredRefinementForNextIteration
               , MPI_CXX_BOOL		 //hasErased
               , MPI_CXX_BOOL		 //hasTriggeredEraseForNextIteration
               , MPI_CXX_BOOL		 //hasChangedVertexOrCellState
               , MPI_CXX_BOOL		 //hasModifiedGridInPreviousIteration
               , MPI_CXX_BOOL		 //isTraversalInverted
               #ifndef MPI2
               , MPI_UB
               #endif
               
            };
            
            int blocklen[Attributes] = {
                 1		 //maxRefinementLevelAllowed
               , 1		 //mergeMode
               , 1		 //sendMode
               , 1		 //_algorithmSection
               , 1		 //hasRefined
               , 1		 //hasTriggeredRefinementForNextIteration
               , 1		 //hasErased
               , 1		 //hasTriggeredEraseForNextIteration
               , 1		 //hasChangedVertexOrCellState
               , 1		 //hasModifiedGridInPreviousIteration
               , 1		 //isTraversalInverted
               #ifndef MPI2
               , 1
               #endif
               
            };
            
            MPI_Aint  disp[Attributes];
            MPI_Aint  base;
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[4] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[4] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[5] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[5] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[6] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[6] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[7] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[7] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[8] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[8] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[9] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[9] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[10] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[10] );
            #endif
            #ifdef MPI2
            for (int i=1; i<Attributes; i++) {
            #else
            for (int i=1; i<Attributes-1; i++) {
            #endif
               assertion1( disp[i] > disp[i-1], i );
            }
            #ifdef MPI2
            for (int i=0; i<Attributes; i++) {
            #else
            for (int i=0; i<Attributes-1; i++) {
            #endif
               disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
               assertion4(disp[i]<static_cast<int>(sizeof(State)), i, disp[i], Attributes, sizeof(State));
            }
            #ifndef MPI2
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]))), 		&disp[11] );
            disp[11] -= base;
            disp[11] += disp[0];
            #endif
            #ifdef MPI2
            MPI_Datatype tmpType; 
            MPI_Aint lowerBound, typeExtent; 
            MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
            MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
            MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &State::Datatype );
            MPI_Type_commit( &State::Datatype );
            #else
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::Datatype);
            MPI_Type_commit( &State::Datatype );
            #endif
            
         }
         {
            State dummyState[2];
            
            #ifdef MPI2
            const int Attributes = 11;
            #else
            const int Attributes = 12;
            #endif
            MPI_Datatype subtypes[Attributes] = {
                 MPI_INT		 //maxRefinementLevelAllowed
               , MPI_INT		 //mergeMode
               , MPI_INT		 //sendMode
               , MPI_INT		 //_algorithmSection
               , MPI_CXX_BOOL		 //hasRefined
               , MPI_CXX_BOOL		 //hasTriggeredRefinementForNextIteration
               , MPI_CXX_BOOL		 //hasErased
               , MPI_CXX_BOOL		 //hasTriggeredEraseForNextIteration
               , MPI_CXX_BOOL		 //hasChangedVertexOrCellState
               , MPI_CXX_BOOL		 //hasModifiedGridInPreviousIteration
               , MPI_CXX_BOOL		 //isTraversalInverted
               #ifndef MPI2
               , MPI_UB
               #endif
               
            };
            
            int blocklen[Attributes] = {
                 1		 //maxRefinementLevelAllowed
               , 1		 //mergeMode
               , 1		 //sendMode
               , 1		 //_algorithmSection
               , 1		 //hasRefined
               , 1		 //hasTriggeredRefinementForNextIteration
               , 1		 //hasErased
               , 1		 //hasTriggeredEraseForNextIteration
               , 1		 //hasChangedVertexOrCellState
               , 1		 //hasModifiedGridInPreviousIteration
               , 1		 //isTraversalInverted
               #ifndef MPI2
               , 1
               #endif
               
            };
            
            MPI_Aint  disp[Attributes];
            MPI_Aint  base;
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[4] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[4] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[5] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[5] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[6] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[6] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[7] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[7] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[8] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[8] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[9] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[9] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[10] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[10] );
            #endif
            #ifdef MPI2
            for (int i=1; i<Attributes; i++) {
            #else
            for (int i=1; i<Attributes-1; i++) {
            #endif
               assertion1( disp[i] > disp[i-1], i );
            }
            #ifdef MPI2
            for (int i=0; i<Attributes; i++) {
            #else
            for (int i=0; i<Attributes-1; i++) {
            #endif
               disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
               assertion4(disp[i]<static_cast<int>(sizeof(State)), i, disp[i], Attributes, sizeof(State));
            }
            #ifndef MPI2
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]))), 		&disp[11] );
            disp[11] -= base;
            disp[11] += disp[0];
            #endif
            #ifdef MPI2
            MPI_Datatype tmpType; 
            MPI_Aint lowerBound, typeExtent; 
            MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
            MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
            MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &State::FullDatatype );
            MPI_Type_commit( &State::FullDatatype );
            #else
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::FullDatatype);
            MPI_Type_commit( &State::FullDatatype );
            #endif
            
         }
         
      }
      
      
      void exahype::records::State::shutdownDatatype() {
         MPI_Type_free( &State::Datatype );
         MPI_Type_free( &State::FullDatatype );
         
      }
      
      void exahype::records::State::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
         _senderDestinationRank = destination;
         
         if (communicateSleep<0) {
         
            const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
            if  (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "was not able to send message exahype::records::State "
               << toString()
               << " to node " << destination
               << ": " << tarch::parallel::MPIReturnValueToString(result);
               _log.error( "send(int)",msg.str() );
            }
            
         }
         else {
         
            MPI_Request* sendRequestHandle = new MPI_Request();
            MPI_Status   status;
            int          flag = 0;
            int          result;
            
            clock_t      timeOutWarning   = -1;
            clock_t      timeOutShutdown  = -1;
            bool         triggeredTimeoutWarning = false;
            
            if (exchangeOnlyAttributesMarkedWithParallelise) {
               result = MPI_Isend(
                  this, 1, Datatype, destination,
                  tag, tarch::parallel::Node::getInstance().getCommunicator(),
                  sendRequestHandle
               );
               
            }
            else {
               result = MPI_Isend(
                  this, 1, FullDatatype, destination,
                  tag, tarch::parallel::Node::getInstance().getCommunicator(),
                  sendRequestHandle
               );
               
            }
            if  (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "was not able to send message exahype::records::State "
               << toString()
               << " to node " << destination
               << ": " << tarch::parallel::MPIReturnValueToString(result);
               _log.error( "send(int)",msg.str() );
            }
            result = MPI_Test( sendRequestHandle, &flag, &status );
            while (!flag) {
               if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
               if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
               result = MPI_Test( sendRequestHandle, &flag, &status );
               if (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "testing for finished send task for exahype::records::State "
                  << toString()
                  << " sent to node " << destination
                  << " failed: " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error("send(int)", msg.str() );
               }
               
               // deadlock aspect
               if (
                  tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                  (clock()>timeOutWarning) &&
                  (!triggeredTimeoutWarning)
               ) {
                  tarch::parallel::Node::getInstance().writeTimeOutWarning(
                  "exahype::records::State",
                  "send(int)", destination,tag,1
                  );
                  triggeredTimeoutWarning = true;
               }
               if (
                  tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                  (clock()>timeOutShutdown)
               ) {
                  tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                  "exahype::records::State",
                  "send(int)", destination,tag,1
                  );
               }
               
            tarch::parallel::Node::getInstance().receiveDanglingMessages();
            usleep(communicateSleep);
            }
            
            delete sendRequestHandle;
            #ifdef Debug
            _log.debug("send(int,int)", "sent " + toString() );
            #endif
            
         }
         
      }
      
      
      
      void exahype::records::State::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
         if (communicateSleep<0) {
         
            MPI_Status  status;
            const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
            _senderDestinationRank = status.MPI_SOURCE;
            if ( result != MPI_SUCCESS ) {
               std::ostringstream msg;
               msg << "failed to start to receive exahype::records::State from node "
               << source << ": " << tarch::parallel::MPIReturnValueToString(result);
               _log.error( "receive(int)", msg.str() );
            }
            
         }
         else {
         
            MPI_Request* sendRequestHandle = new MPI_Request();
            MPI_Status   status;
            int          flag = 0;
            int          result;
            
            clock_t      timeOutWarning   = -1;
            clock_t      timeOutShutdown  = -1;
            bool         triggeredTimeoutWarning = false;
            
            if (exchangeOnlyAttributesMarkedWithParallelise) {
               result = MPI_Irecv(
                  this, 1, Datatype, source, tag,
                  tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
               );
               
            }
            else {
               result = MPI_Irecv(
                  this, 1, FullDatatype, source, tag,
                  tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
               );
               
            }
            if ( result != MPI_SUCCESS ) {
               std::ostringstream msg;
               msg << "failed to start to receive exahype::records::State from node "
               << source << ": " << tarch::parallel::MPIReturnValueToString(result);
               _log.error( "receive(int)", msg.str() );
            }
            
            result = MPI_Test( sendRequestHandle, &flag, &status );
            while (!flag) {
               if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
               if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
               result = MPI_Test( sendRequestHandle, &flag, &status );
               if (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "testing for finished receive task for exahype::records::State failed: "
                  << tarch::parallel::MPIReturnValueToString(result);
                  _log.error("receive(int)", msg.str() );
               }
               
               // deadlock aspect
               if (
                  tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                  (clock()>timeOutWarning) &&
                  (!triggeredTimeoutWarning)
               ) {
                  tarch::parallel::Node::getInstance().writeTimeOutWarning(
                  "exahype::records::State",
                  "receive(int)", source,tag,1
                  );
                  triggeredTimeoutWarning = true;
               }
               if (
                  tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                  (clock()>timeOutShutdown)
               ) {
                  tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                  "exahype::records::State",
                  "receive(int)", source,tag,1
                  );
               }
               tarch::parallel::Node::getInstance().receiveDanglingMessages();
               usleep(communicateSleep);
               
            }
            
            delete sendRequestHandle;
            
            _senderDestinationRank = status.MPI_SOURCE;
            #ifdef Debug
            _log.debug("receive(int,int)", "received " + toString() ); 
            #endif
            
         }
         
      }
      
      
      
      bool exahype::records::State::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
         MPI_Status status;
         int  flag        = 0;
         MPI_Iprobe(
            MPI_ANY_SOURCE, tag,
            tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
         );
         if (flag) {
            int  messageCounter;
            if (exchangeOnlyAttributesMarkedWithParallelise) {
               MPI_Get_count(&status, Datatype, &messageCounter);
            }
            else {
               MPI_Get_count(&status, FullDatatype, &messageCounter);
            }
            return messageCounter > 0;
         }
         else return false;
         
      }
      
      int exahype::records::State::getSenderRank() const {
         assertion( _senderDestinationRank!=-1 );
         return _senderDestinationRank;
         
      }
   #endif
   
   
   exahype::records::StatePacked::PersistentRecords::PersistentRecords() {
      if ((6 >= (8 * sizeof(short int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((6 < (8 * sizeof(short int))));
      
   }
   
   
   exahype::records::StatePacked::PersistentRecords::PersistentRecords(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
   _maxRefinementLevelAllowed(maxRefinementLevelAllowed),
   _mergeMode(mergeMode),
   _sendMode(sendMode),
   _algorithmSection(algorithmSection),
   _isTraversalInverted(isTraversalInverted) {
      setHasRefined(hasRefined);
      setHasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration);
      setHasErased(hasErased);
      setHasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration);
      setHasChangedVertexOrCellState(hasChangedVertexOrCellState);
      setHasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration);
      if ((6 >= (8 * sizeof(short int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((6 < (8 * sizeof(short int))));
      
   }
   
   exahype::records::StatePacked::StatePacked() {
      if ((6 >= (8 * sizeof(short int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((6 < (8 * sizeof(short int))));
      
   }
   
   
   exahype::records::StatePacked::StatePacked(const PersistentRecords& persistentRecords):
   _persistentRecords(persistentRecords._maxRefinementLevelAllowed, persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._algorithmSection, persistentRecords.getHasRefined(), persistentRecords.getHasTriggeredRefinementForNextIteration(), persistentRecords.getHasErased(), persistentRecords.getHasTriggeredEraseForNextIteration(), persistentRecords.getHasChangedVertexOrCellState(), persistentRecords.getHasModifiedGridInPreviousIteration(), persistentRecords._isTraversalInverted) {
      if ((6 >= (8 * sizeof(short int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((6 < (8 * sizeof(short int))));
      
   }
   
   
   exahype::records::StatePacked::StatePacked(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
   _persistentRecords(maxRefinementLevelAllowed, mergeMode, sendMode, algorithmSection, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted) {
      if ((6 >= (8 * sizeof(short int)))) {
         std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
         std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
         std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
      }
      assertion((6 < (8 * sizeof(short int))));
      
   }
   
   
   exahype::records::StatePacked::~StatePacked() { }
   
   std::string exahype::records::StatePacked::toString(const MergeMode& param) {
      return exahype::records::State::toString(param);
   }
   
   std::string exahype::records::StatePacked::getMergeModeMapping() {
      return exahype::records::State::getMergeModeMapping();
   }
   
   std::string exahype::records::StatePacked::toString(const SendMode& param) {
      return exahype::records::State::toString(param);
   }
   
   std::string exahype::records::StatePacked::getSendModeMapping() {
      return exahype::records::State::getSendModeMapping();
   }
   
   std::string exahype::records::StatePacked::toString(const AlgorithmSection& param) {
      return exahype::records::State::toString(param);
   }
   
   std::string exahype::records::StatePacked::getAlgorithmSectionMapping() {
      return exahype::records::State::getAlgorithmSectionMapping();
   }
   
   
   
   std::string exahype::records::StatePacked::toString() const {
      std::ostringstream stringstr;
      toString(stringstr);
      return stringstr.str();
   }
   
   void exahype::records::StatePacked::toString (std::ostream& out) const {
      out << "("; 
      out << "maxRefinementLevelAllowed:" << getMaxRefinementLevelAllowed();
      out << ",";
      out << "mergeMode:" << toString(getMergeMode());
      out << ",";
      out << "sendMode:" << toString(getSendMode());
      out << ",";
      out << "_algorithmSection:" << toString(getAlgorithmSection());
      out << ",";
      out << "hasRefined:" << getHasRefined();
      out << ",";
      out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
      out << ",";
      out << "hasErased:" << getHasErased();
      out << ",";
      out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
      out << ",";
      out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
      out << ",";
      out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
      out << ",";
      out << "isTraversalInverted:" << getIsTraversalInverted();
      out <<  ")";
   }
   
   
   exahype::records::StatePacked::PersistentRecords exahype::records::StatePacked::getPersistentRecords() const {
      return _persistentRecords;
   }
   
   exahype::records::State exahype::records::StatePacked::convert() const{
      return State(
         getMaxRefinementLevelAllowed(),
         getMergeMode(),
         getSendMode(),
         getAlgorithmSection(),
         getHasRefined(),
         getHasTriggeredRefinementForNextIteration(),
         getHasErased(),
         getHasTriggeredEraseForNextIteration(),
         getHasChangedVertexOrCellState(),
         getHasModifiedGridInPreviousIteration(),
         getIsTraversalInverted()
      );
   }
   
   #ifdef Parallel
      tarch::logging::Log exahype::records::StatePacked::_log( "exahype::records::StatePacked" );
      
      MPI_Datatype exahype::records::StatePacked::Datatype = 0;
      MPI_Datatype exahype::records::StatePacked::FullDatatype = 0;
      
      
      void exahype::records::StatePacked::initDatatype() {
         {
            StatePacked dummyStatePacked[2];
            
            #ifdef MPI2
            const int Attributes = 6;
            #else
            const int Attributes = 7;
            #endif
            MPI_Datatype subtypes[Attributes] = {
                 MPI_INT		 //maxRefinementLevelAllowed
               , MPI_INT		 //mergeMode
               , MPI_INT		 //sendMode
               , MPI_INT		 //_algorithmSection
               , MPI_CXX_BOOL		 //isTraversalInverted
               , MPI_SHORT		 //_packedRecords0
               #ifndef MPI2
               , MPI_UB
               #endif
               
            };
            
            int blocklen[Attributes] = {
                 1		 //maxRefinementLevelAllowed
               , 1		 //mergeMode
               , 1		 //sendMode
               , 1		 //_algorithmSection
               , 1		 //isTraversalInverted
               , 1		 //_packedRecords0
               #ifndef MPI2
               , 1
               #endif
               
            };
            
            MPI_Aint  disp[Attributes];
            MPI_Aint  base;
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[4] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[4] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[5] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[5] );
            #endif
            #ifdef MPI2
            for (int i=1; i<Attributes; i++) {
            #else
            for (int i=1; i<Attributes-1; i++) {
            #endif
               assertion1( disp[i] > disp[i-1], i );
            }
            #ifdef MPI2
            for (int i=0; i<Attributes; i++) {
            #else
            for (int i=0; i<Attributes-1; i++) {
            #endif
               disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
               assertion4(disp[i]<static_cast<int>(sizeof(StatePacked)), i, disp[i], Attributes, sizeof(StatePacked));
            }
            #ifndef MPI2
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]))), 		&disp[6] );
            disp[6] -= base;
            disp[6] += disp[0];
            #endif
            #ifdef MPI2
            MPI_Datatype tmpType; 
            MPI_Aint lowerBound, typeExtent; 
            MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
            MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
            MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &StatePacked::Datatype );
            MPI_Type_commit( &StatePacked::Datatype );
            #else
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::Datatype);
            MPI_Type_commit( &StatePacked::Datatype );
            #endif
            
         }
         {
            StatePacked dummyStatePacked[2];
            
            #ifdef MPI2
            const int Attributes = 6;
            #else
            const int Attributes = 7;
            #endif
            MPI_Datatype subtypes[Attributes] = {
                 MPI_INT		 //maxRefinementLevelAllowed
               , MPI_INT		 //mergeMode
               , MPI_INT		 //sendMode
               , MPI_INT		 //_algorithmSection
               , MPI_CXX_BOOL		 //isTraversalInverted
               , MPI_SHORT		 //_packedRecords0
               #ifndef MPI2
               , MPI_UB
               #endif
               
            };
            
            int blocklen[Attributes] = {
                 1		 //maxRefinementLevelAllowed
               , 1		 //mergeMode
               , 1		 //sendMode
               , 1		 //_algorithmSection
               , 1		 //isTraversalInverted
               , 1		 //_packedRecords0
               #ifndef MPI2
               , 1
               #endif
               
            };
            
            MPI_Aint  disp[Attributes];
            MPI_Aint  base;
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[4] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[4] );
            #endif
            #ifdef MPI2
            MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[5] );
            #else
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[5] );
            #endif
            #ifdef MPI2
            for (int i=1; i<Attributes; i++) {
            #else
            for (int i=1; i<Attributes-1; i++) {
            #endif
               assertion1( disp[i] > disp[i-1], i );
            }
            #ifdef MPI2
            for (int i=0; i<Attributes; i++) {
            #else
            for (int i=0; i<Attributes-1; i++) {
            #endif
               disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
               assertion4(disp[i]<static_cast<int>(sizeof(StatePacked)), i, disp[i], Attributes, sizeof(StatePacked));
            }
            #ifndef MPI2
            MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]))), 		&disp[6] );
            disp[6] -= base;
            disp[6] += disp[0];
            #endif
            #ifdef MPI2
            MPI_Datatype tmpType; 
            MPI_Aint lowerBound, typeExtent; 
            MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
            MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
            MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &StatePacked::FullDatatype );
            MPI_Type_commit( &StatePacked::FullDatatype );
            #else
            MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::FullDatatype);
            MPI_Type_commit( &StatePacked::FullDatatype );
            #endif
            
         }
         
      }
      
      
      void exahype::records::StatePacked::shutdownDatatype() {
         MPI_Type_free( &StatePacked::Datatype );
         MPI_Type_free( &StatePacked::FullDatatype );
         
      }
      
      void exahype::records::StatePacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
         _senderDestinationRank = destination;
         
         if (communicateSleep<0) {
         
            const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
            if  (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "was not able to send message exahype::records::StatePacked "
               << toString()
               << " to node " << destination
               << ": " << tarch::parallel::MPIReturnValueToString(result);
               _log.error( "send(int)",msg.str() );
            }
            
         }
         else {
         
            MPI_Request* sendRequestHandle = new MPI_Request();
            MPI_Status   status;
            int          flag = 0;
            int          result;
            
            clock_t      timeOutWarning   = -1;
            clock_t      timeOutShutdown  = -1;
            bool         triggeredTimeoutWarning = false;
            
            if (exchangeOnlyAttributesMarkedWithParallelise) {
               result = MPI_Isend(
                  this, 1, Datatype, destination,
                  tag, tarch::parallel::Node::getInstance().getCommunicator(),
                  sendRequestHandle
               );
               
            }
            else {
               result = MPI_Isend(
                  this, 1, FullDatatype, destination,
                  tag, tarch::parallel::Node::getInstance().getCommunicator(),
                  sendRequestHandle
               );
               
            }
            if  (result!=MPI_SUCCESS) {
               std::ostringstream msg;
               msg << "was not able to send message exahype::records::StatePacked "
               << toString()
               << " to node " << destination
               << ": " << tarch::parallel::MPIReturnValueToString(result);
               _log.error( "send(int)",msg.str() );
            }
            result = MPI_Test( sendRequestHandle, &flag, &status );
            while (!flag) {
               if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
               if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
               result = MPI_Test( sendRequestHandle, &flag, &status );
               if (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "testing for finished send task for exahype::records::StatePacked "
                  << toString()
                  << " sent to node " << destination
                  << " failed: " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error("send(int)", msg.str() );
               }
               
               // deadlock aspect
               if (
                  tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                  (clock()>timeOutWarning) &&
                  (!triggeredTimeoutWarning)
               ) {
                  tarch::parallel::Node::getInstance().writeTimeOutWarning(
                  "exahype::records::StatePacked",
                  "send(int)", destination,tag,1
                  );
                  triggeredTimeoutWarning = true;
               }
               if (
                  tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                  (clock()>timeOutShutdown)
               ) {
                  tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                  "exahype::records::StatePacked",
                  "send(int)", destination,tag,1
                  );
               }
               
            tarch::parallel::Node::getInstance().receiveDanglingMessages();
            usleep(communicateSleep);
            }
            
            delete sendRequestHandle;
            #ifdef Debug
            _log.debug("send(int,int)", "sent " + toString() );
            #endif
            
         }
         
      }
      
      
      
      void exahype::records::StatePacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
         if (communicateSleep<0) {
         
            MPI_Status  status;
            const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
            _senderDestinationRank = status.MPI_SOURCE;
            if ( result != MPI_SUCCESS ) {
               std::ostringstream msg;
               msg << "failed to start to receive exahype::records::StatePacked from node "
               << source << ": " << tarch::parallel::MPIReturnValueToString(result);
               _log.error( "receive(int)", msg.str() );
            }
            
         }
         else {
         
            MPI_Request* sendRequestHandle = new MPI_Request();
            MPI_Status   status;
            int          flag = 0;
            int          result;
            
            clock_t      timeOutWarning   = -1;
            clock_t      timeOutShutdown  = -1;
            bool         triggeredTimeoutWarning = false;
            
            if (exchangeOnlyAttributesMarkedWithParallelise) {
               result = MPI_Irecv(
                  this, 1, Datatype, source, tag,
                  tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
               );
               
            }
            else {
               result = MPI_Irecv(
                  this, 1, FullDatatype, source, tag,
                  tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
               );
               
            }
            if ( result != MPI_SUCCESS ) {
               std::ostringstream msg;
               msg << "failed to start to receive exahype::records::StatePacked from node "
               << source << ": " << tarch::parallel::MPIReturnValueToString(result);
               _log.error( "receive(int)", msg.str() );
            }
            
            result = MPI_Test( sendRequestHandle, &flag, &status );
            while (!flag) {
               if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
               if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
               result = MPI_Test( sendRequestHandle, &flag, &status );
               if (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "testing for finished receive task for exahype::records::StatePacked failed: "
                  << tarch::parallel::MPIReturnValueToString(result);
                  _log.error("receive(int)", msg.str() );
               }
               
               // deadlock aspect
               if (
                  tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                  (clock()>timeOutWarning) &&
                  (!triggeredTimeoutWarning)
               ) {
                  tarch::parallel::Node::getInstance().writeTimeOutWarning(
                  "exahype::records::StatePacked",
                  "receive(int)", source,tag,1
                  );
                  triggeredTimeoutWarning = true;
               }
               if (
                  tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                  (clock()>timeOutShutdown)
               ) {
                  tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                  "exahype::records::StatePacked",
                  "receive(int)", source,tag,1
                  );
               }
               tarch::parallel::Node::getInstance().receiveDanglingMessages();
               usleep(communicateSleep);
               
            }
            
            delete sendRequestHandle;
            
            _senderDestinationRank = status.MPI_SOURCE;
            #ifdef Debug
            _log.debug("receive(int,int)", "received " + toString() ); 
            #endif
            
         }
         
      }
      
      
      
      bool exahype::records::StatePacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
         MPI_Status status;
         int  flag        = 0;
         MPI_Iprobe(
            MPI_ANY_SOURCE, tag,
            tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
         );
         if (flag) {
            int  messageCounter;
            if (exchangeOnlyAttributesMarkedWithParallelise) {
               MPI_Get_count(&status, Datatype, &messageCounter);
            }
            else {
               MPI_Get_count(&status, FullDatatype, &messageCounter);
            }
            return messageCounter > 0;
         }
         else return false;
         
      }
      
      int exahype::records::StatePacked::getSenderRank() const {
         assertion( _senderDestinationRank!=-1 );
         return _senderDestinationRank;
         
      }
   #endif
   
   
   #elif defined(TrackGridStatistics) && defined(Parallel)
      exahype::records::State::PersistentRecords::PersistentRecords() {
         
      }
      
      
      exahype::records::State::PersistentRecords::PersistentRecords(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
      _maxRefinementLevelAllowed(maxRefinementLevelAllowed),
      _mergeMode(mergeMode),
      _sendMode(sendMode),
      _algorithmSection(algorithmSection),
      _minMeshWidth(minMeshWidth),
      _maxMeshWidth(maxMeshWidth),
      _numberOfInnerVertices(numberOfInnerVertices),
      _numberOfBoundaryVertices(numberOfBoundaryVertices),
      _numberOfOuterVertices(numberOfOuterVertices),
      _numberOfInnerCells(numberOfInnerCells),
      _numberOfOuterCells(numberOfOuterCells),
      _numberOfInnerLeafVertices(numberOfInnerLeafVertices),
      _numberOfBoundaryLeafVertices(numberOfBoundaryLeafVertices),
      _numberOfOuterLeafVertices(numberOfOuterLeafVertices),
      _numberOfInnerLeafCells(numberOfInnerLeafCells),
      _numberOfOuterLeafCells(numberOfOuterLeafCells),
      _maxLevel(maxLevel),
      _hasRefined(hasRefined),
      _hasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration),
      _hasErased(hasErased),
      _hasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration),
      _hasChangedVertexOrCellState(hasChangedVertexOrCellState),
      _hasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration),
      _isTraversalInverted(isTraversalInverted),
      _reduceStateAndCell(reduceStateAndCell),
      _couldNotEraseDueToDecompositionFlag(couldNotEraseDueToDecompositionFlag),
      _subWorkerIsInvolvedInJoinOrFork(subWorkerIsInvolvedInJoinOrFork) {
         
      }
      
      exahype::records::State::State() {
         
      }
      
      
      exahype::records::State::State(const PersistentRecords& persistentRecords):
      _persistentRecords(persistentRecords._maxRefinementLevelAllowed, persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._algorithmSection, persistentRecords._minMeshWidth, persistentRecords._maxMeshWidth, persistentRecords._numberOfInnerVertices, persistentRecords._numberOfBoundaryVertices, persistentRecords._numberOfOuterVertices, persistentRecords._numberOfInnerCells, persistentRecords._numberOfOuterCells, persistentRecords._numberOfInnerLeafVertices, persistentRecords._numberOfBoundaryLeafVertices, persistentRecords._numberOfOuterLeafVertices, persistentRecords._numberOfInnerLeafCells, persistentRecords._numberOfOuterLeafCells, persistentRecords._maxLevel, persistentRecords._hasRefined, persistentRecords._hasTriggeredRefinementForNextIteration, persistentRecords._hasErased, persistentRecords._hasTriggeredEraseForNextIteration, persistentRecords._hasChangedVertexOrCellState, persistentRecords._hasModifiedGridInPreviousIteration, persistentRecords._isTraversalInverted, persistentRecords._reduceStateAndCell, persistentRecords._couldNotEraseDueToDecompositionFlag, persistentRecords._subWorkerIsInvolvedInJoinOrFork) {
         
      }
      
      
      exahype::records::State::State(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
      _persistentRecords(maxRefinementLevelAllowed, mergeMode, sendMode, algorithmSection, minMeshWidth, maxMeshWidth, numberOfInnerVertices, numberOfBoundaryVertices, numberOfOuterVertices, numberOfInnerCells, numberOfOuterCells, numberOfInnerLeafVertices, numberOfBoundaryLeafVertices, numberOfOuterLeafVertices, numberOfInnerLeafCells, numberOfOuterLeafCells, maxLevel, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted, reduceStateAndCell, couldNotEraseDueToDecompositionFlag, subWorkerIsInvolvedInJoinOrFork) {
         
      }
      
      
      exahype::records::State::~State() { }
      
      std::string exahype::records::State::toString(const AlgorithmSection& param) {
         switch (param) {
            case TimeStepping: return "TimeStepping";
            case LimiterStatusSpreading: return "LimiterStatusSpreading";
            case MeshRefinement: return "MeshRefinement";
            case MeshRefinementOrLocalOrGlobalRecomputation: return "MeshRefinementOrLocalOrGlobalRecomputation";
            case LocalRecomputationAllSend: return "LocalRecomputationAllSend";
            case MeshRefinementOrGlobalRecomputation: return "MeshRefinementOrGlobalRecomputation";
            case MeshRefinementOrGlobalRecomputationAllSend: return "MeshRefinementOrGlobalRecomputationAllSend";
            case PredictionRerunAllSend: return "PredictionRerunAllSend";
         }
         return "undefined";
      }
      
      std::string exahype::records::State::getAlgorithmSectionMapping() {
         return "AlgorithmSection(TimeStepping=0,LimiterStatusSpreading=1,MeshRefinement=2,MeshRefinementOrLocalOrGlobalRecomputation=3,LocalRecomputationAllSend=4,MeshRefinementOrGlobalRecomputation=5,MeshRefinementOrGlobalRecomputationAllSend=6,PredictionRerunAllSend=7)";
      }
      std::string exahype::records::State::toString(const MergeMode& param) {
         switch (param) {
            case MergeNothing: return "MergeNothing";
            case BroadcastAndMergeTimeStepData: return "BroadcastAndMergeTimeStepData";
            case MergeFaceData: return "MergeFaceData";
            case DropFaceData: return "DropFaceData";
            case BroadcastAndMergeTimeStepDataAndMergeFaceData: return "BroadcastAndMergeTimeStepDataAndMergeFaceData";
            case BroadcastAndMergeTimeStepDataAndDropFaceData: return "BroadcastAndMergeTimeStepDataAndDropFaceData";
         }
         return "undefined";
      }
      
      std::string exahype::records::State::getMergeModeMapping() {
         return "MergeMode(MergeNothing=0,BroadcastAndMergeTimeStepData=1,MergeFaceData=2,DropFaceData=3,BroadcastAndMergeTimeStepDataAndMergeFaceData=4,BroadcastAndMergeTimeStepDataAndDropFaceData=5)";
      }
      std::string exahype::records::State::toString(const SendMode& param) {
         switch (param) {
            case SendNothing: return "SendNothing";
            case ReduceAndMergeTimeStepData: return "ReduceAndMergeTimeStepData";
            case SendFaceData: return "SendFaceData";
            case ReduceAndMergeTimeStepDataAndSendFaceData: return "ReduceAndMergeTimeStepDataAndSendFaceData";
         }
         return "undefined";
      }
      
      std::string exahype::records::State::getSendModeMapping() {
         return "SendMode(SendNothing=0,ReduceAndMergeTimeStepData=1,SendFaceData=2,ReduceAndMergeTimeStepDataAndSendFaceData=3)";
      }
      
      
      std::string exahype::records::State::toString() const {
         std::ostringstream stringstr;
         toString(stringstr);
         return stringstr.str();
      }
      
      void exahype::records::State::toString (std::ostream& out) const {
         out << "("; 
         out << "maxRefinementLevelAllowed:" << getMaxRefinementLevelAllowed();
         out << ",";
         out << "mergeMode:" << toString(getMergeMode());
         out << ",";
         out << "sendMode:" << toString(getSendMode());
         out << ",";
         out << "_algorithmSection:" << toString(getAlgorithmSection());
         out << ",";
         out << "minMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMinMeshWidth(i) << ",";
   }
   out << getMinMeshWidth(DIMENSIONS-1) << "]";
         out << ",";
         out << "maxMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMaxMeshWidth(i) << ",";
   }
   out << getMaxMeshWidth(DIMENSIONS-1) << "]";
         out << ",";
         out << "numberOfInnerVertices:" << getNumberOfInnerVertices();
         out << ",";
         out << "numberOfBoundaryVertices:" << getNumberOfBoundaryVertices();
         out << ",";
         out << "numberOfOuterVertices:" << getNumberOfOuterVertices();
         out << ",";
         out << "numberOfInnerCells:" << getNumberOfInnerCells();
         out << ",";
         out << "numberOfOuterCells:" << getNumberOfOuterCells();
         out << ",";
         out << "numberOfInnerLeafVertices:" << getNumberOfInnerLeafVertices();
         out << ",";
         out << "numberOfBoundaryLeafVertices:" << getNumberOfBoundaryLeafVertices();
         out << ",";
         out << "numberOfOuterLeafVertices:" << getNumberOfOuterLeafVertices();
         out << ",";
         out << "numberOfInnerLeafCells:" << getNumberOfInnerLeafCells();
         out << ",";
         out << "numberOfOuterLeafCells:" << getNumberOfOuterLeafCells();
         out << ",";
         out << "maxLevel:" << getMaxLevel();
         out << ",";
         out << "hasRefined:" << getHasRefined();
         out << ",";
         out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
         out << ",";
         out << "hasErased:" << getHasErased();
         out << ",";
         out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
         out << ",";
         out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
         out << ",";
         out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
         out << ",";
         out << "isTraversalInverted:" << getIsTraversalInverted();
         out << ",";
         out << "reduceStateAndCell:" << getReduceStateAndCell();
         out << ",";
         out << "couldNotEraseDueToDecompositionFlag:" << getCouldNotEraseDueToDecompositionFlag();
         out << ",";
         out << "subWorkerIsInvolvedInJoinOrFork:" << getSubWorkerIsInvolvedInJoinOrFork();
         out <<  ")";
      }
      
      
      exahype::records::State::PersistentRecords exahype::records::State::getPersistentRecords() const {
         return _persistentRecords;
      }
      
      exahype::records::StatePacked exahype::records::State::convert() const{
         return StatePacked(
            getMaxRefinementLevelAllowed(),
            getMergeMode(),
            getSendMode(),
            getAlgorithmSection(),
            getMinMeshWidth(),
            getMaxMeshWidth(),
            getNumberOfInnerVertices(),
            getNumberOfBoundaryVertices(),
            getNumberOfOuterVertices(),
            getNumberOfInnerCells(),
            getNumberOfOuterCells(),
            getNumberOfInnerLeafVertices(),
            getNumberOfBoundaryLeafVertices(),
            getNumberOfOuterLeafVertices(),
            getNumberOfInnerLeafCells(),
            getNumberOfOuterLeafCells(),
            getMaxLevel(),
            getHasRefined(),
            getHasTriggeredRefinementForNextIteration(),
            getHasErased(),
            getHasTriggeredEraseForNextIteration(),
            getHasChangedVertexOrCellState(),
            getHasModifiedGridInPreviousIteration(),
            getIsTraversalInverted(),
            getReduceStateAndCell(),
            getCouldNotEraseDueToDecompositionFlag(),
            getSubWorkerIsInvolvedInJoinOrFork()
         );
      }
      
      #ifdef Parallel
         tarch::logging::Log exahype::records::State::_log( "exahype::records::State" );
         
         MPI_Datatype exahype::records::State::Datatype = 0;
         MPI_Datatype exahype::records::State::FullDatatype = 0;
         
         
         void exahype::records::State::initDatatype() {
            {
               State dummyState[2];
               
               #ifdef MPI2
               const int Attributes = 27;
               #else
               const int Attributes = 28;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //maxRefinementLevelAllowed
                  , MPI_INT		 //mergeMode
                  , MPI_INT		 //sendMode
                  , MPI_INT		 //_algorithmSection
                  , MPI_DOUBLE		 //minMeshWidth
                  , MPI_DOUBLE		 //maxMeshWidth
                  , MPI_DOUBLE		 //numberOfInnerVertices
                  , MPI_DOUBLE		 //numberOfBoundaryVertices
                  , MPI_DOUBLE		 //numberOfOuterVertices
                  , MPI_DOUBLE		 //numberOfInnerCells
                  , MPI_DOUBLE		 //numberOfOuterCells
                  , MPI_DOUBLE		 //numberOfInnerLeafVertices
                  , MPI_DOUBLE		 //numberOfBoundaryLeafVertices
                  , MPI_DOUBLE		 //numberOfOuterLeafVertices
                  , MPI_DOUBLE		 //numberOfInnerLeafCells
                  , MPI_DOUBLE		 //numberOfOuterLeafCells
                  , MPI_INT		 //maxLevel
                  , MPI_CXX_BOOL		 //hasRefined
                  , MPI_CXX_BOOL		 //hasTriggeredRefinementForNextIteration
                  , MPI_CXX_BOOL		 //hasErased
                  , MPI_CXX_BOOL		 //hasTriggeredEraseForNextIteration
                  , MPI_CXX_BOOL		 //hasChangedVertexOrCellState
                  , MPI_CXX_BOOL		 //hasModifiedGridInPreviousIteration
                  , MPI_CXX_BOOL		 //isTraversalInverted
                  , MPI_CXX_BOOL		 //reduceStateAndCell
                  , MPI_CXX_BOOL		 //couldNotEraseDueToDecompositionFlag
                  , MPI_CXX_BOOL		 //subWorkerIsInvolvedInJoinOrFork
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //maxRefinementLevelAllowed
                  , 1		 //mergeMode
                  , 1		 //sendMode
                  , 1		 //_algorithmSection
                  , DIMENSIONS		 //minMeshWidth
                  , DIMENSIONS		 //maxMeshWidth
                  , 1		 //numberOfInnerVertices
                  , 1		 //numberOfBoundaryVertices
                  , 1		 //numberOfOuterVertices
                  , 1		 //numberOfInnerCells
                  , 1		 //numberOfOuterCells
                  , 1		 //numberOfInnerLeafVertices
                  , 1		 //numberOfBoundaryLeafVertices
                  , 1		 //numberOfOuterLeafVertices
                  , 1		 //numberOfInnerLeafCells
                  , 1		 //numberOfOuterLeafCells
                  , 1		 //maxLevel
                  , 1		 //hasRefined
                  , 1		 //hasTriggeredRefinementForNextIteration
                  , 1		 //hasErased
                  , 1		 //hasTriggeredEraseForNextIteration
                  , 1		 //hasChangedVertexOrCellState
                  , 1		 //hasModifiedGridInPreviousIteration
                  , 1		 //isTraversalInverted
                  , 1		 //reduceStateAndCell
                  , 1		 //couldNotEraseDueToDecompositionFlag
                  , 1		 //subWorkerIsInvolvedInJoinOrFork
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[18] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[19] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[19] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[20] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[20] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[21] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[21] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[22] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[22] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[23] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[23] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reduceStateAndCell))), 		&disp[24] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reduceStateAndCell))), 		&disp[24] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._couldNotEraseDueToDecompositionFlag))), 		&disp[25] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._couldNotEraseDueToDecompositionFlag))), 		&disp[25] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._subWorkerIsInvolvedInJoinOrFork))), 		&disp[26] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._subWorkerIsInvolvedInJoinOrFork))), 		&disp[26] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(State)), i, disp[i], Attributes, sizeof(State));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]))), 		&disp[27] );
               disp[27] -= base;
               disp[27] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &State::Datatype );
               MPI_Type_commit( &State::Datatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::Datatype);
               MPI_Type_commit( &State::Datatype );
               #endif
               
            }
            {
               State dummyState[2];
               
               #ifdef MPI2
               const int Attributes = 27;
               #else
               const int Attributes = 28;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //maxRefinementLevelAllowed
                  , MPI_INT		 //mergeMode
                  , MPI_INT		 //sendMode
                  , MPI_INT		 //_algorithmSection
                  , MPI_DOUBLE		 //minMeshWidth
                  , MPI_DOUBLE		 //maxMeshWidth
                  , MPI_DOUBLE		 //numberOfInnerVertices
                  , MPI_DOUBLE		 //numberOfBoundaryVertices
                  , MPI_DOUBLE		 //numberOfOuterVertices
                  , MPI_DOUBLE		 //numberOfInnerCells
                  , MPI_DOUBLE		 //numberOfOuterCells
                  , MPI_DOUBLE		 //numberOfInnerLeafVertices
                  , MPI_DOUBLE		 //numberOfBoundaryLeafVertices
                  , MPI_DOUBLE		 //numberOfOuterLeafVertices
                  , MPI_DOUBLE		 //numberOfInnerLeafCells
                  , MPI_DOUBLE		 //numberOfOuterLeafCells
                  , MPI_INT		 //maxLevel
                  , MPI_CXX_BOOL		 //hasRefined
                  , MPI_CXX_BOOL		 //hasTriggeredRefinementForNextIteration
                  , MPI_CXX_BOOL		 //hasErased
                  , MPI_CXX_BOOL		 //hasTriggeredEraseForNextIteration
                  , MPI_CXX_BOOL		 //hasChangedVertexOrCellState
                  , MPI_CXX_BOOL		 //hasModifiedGridInPreviousIteration
                  , MPI_CXX_BOOL		 //isTraversalInverted
                  , MPI_CXX_BOOL		 //reduceStateAndCell
                  , MPI_CXX_BOOL		 //couldNotEraseDueToDecompositionFlag
                  , MPI_CXX_BOOL		 //subWorkerIsInvolvedInJoinOrFork
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //maxRefinementLevelAllowed
                  , 1		 //mergeMode
                  , 1		 //sendMode
                  , 1		 //_algorithmSection
                  , DIMENSIONS		 //minMeshWidth
                  , DIMENSIONS		 //maxMeshWidth
                  , 1		 //numberOfInnerVertices
                  , 1		 //numberOfBoundaryVertices
                  , 1		 //numberOfOuterVertices
                  , 1		 //numberOfInnerCells
                  , 1		 //numberOfOuterCells
                  , 1		 //numberOfInnerLeafVertices
                  , 1		 //numberOfBoundaryLeafVertices
                  , 1		 //numberOfOuterLeafVertices
                  , 1		 //numberOfInnerLeafCells
                  , 1		 //numberOfOuterLeafCells
                  , 1		 //maxLevel
                  , 1		 //hasRefined
                  , 1		 //hasTriggeredRefinementForNextIteration
                  , 1		 //hasErased
                  , 1		 //hasTriggeredEraseForNextIteration
                  , 1		 //hasChangedVertexOrCellState
                  , 1		 //hasModifiedGridInPreviousIteration
                  , 1		 //isTraversalInverted
                  , 1		 //reduceStateAndCell
                  , 1		 //couldNotEraseDueToDecompositionFlag
                  , 1		 //subWorkerIsInvolvedInJoinOrFork
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[18] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[19] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[19] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[20] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[20] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[21] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[21] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[22] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[22] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[23] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[23] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reduceStateAndCell))), 		&disp[24] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reduceStateAndCell))), 		&disp[24] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._couldNotEraseDueToDecompositionFlag))), 		&disp[25] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._couldNotEraseDueToDecompositionFlag))), 		&disp[25] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._subWorkerIsInvolvedInJoinOrFork))), 		&disp[26] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._subWorkerIsInvolvedInJoinOrFork))), 		&disp[26] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(State)), i, disp[i], Attributes, sizeof(State));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]))), 		&disp[27] );
               disp[27] -= base;
               disp[27] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &State::FullDatatype );
               MPI_Type_commit( &State::FullDatatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::FullDatatype);
               MPI_Type_commit( &State::FullDatatype );
               #endif
               
            }
            
         }
         
         
         void exahype::records::State::shutdownDatatype() {
            MPI_Type_free( &State::Datatype );
            MPI_Type_free( &State::FullDatatype );
            
         }
         
         void exahype::records::State::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
            _senderDestinationRank = destination;
            
            if (communicateSleep<0) {
            
               const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
               if  (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "was not able to send message exahype::records::State "
                  << toString()
                  << " to node " << destination
                  << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "send(int)",msg.str() );
               }
               
            }
            else {
            
               MPI_Request* sendRequestHandle = new MPI_Request();
               MPI_Status   status;
               int          flag = 0;
               int          result;
               
               clock_t      timeOutWarning   = -1;
               clock_t      timeOutShutdown  = -1;
               bool         triggeredTimeoutWarning = false;
               
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  result = MPI_Isend(
                     this, 1, Datatype, destination,
                     tag, tarch::parallel::Node::getInstance().getCommunicator(),
                     sendRequestHandle
                  );
                  
               }
               else {
                  result = MPI_Isend(
                     this, 1, FullDatatype, destination,
                     tag, tarch::parallel::Node::getInstance().getCommunicator(),
                     sendRequestHandle
                  );
                  
               }
               if  (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "was not able to send message exahype::records::State "
                  << toString()
                  << " to node " << destination
                  << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "send(int)",msg.str() );
               }
               result = MPI_Test( sendRequestHandle, &flag, &status );
               while (!flag) {
                  if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
                  if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
                  result = MPI_Test( sendRequestHandle, &flag, &status );
                  if (result!=MPI_SUCCESS) {
                     std::ostringstream msg;
                     msg << "testing for finished send task for exahype::records::State "
                     << toString()
                     << " sent to node " << destination
                     << " failed: " << tarch::parallel::MPIReturnValueToString(result);
                     _log.error("send(int)", msg.str() );
                  }
                  
                  // deadlock aspect
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                     (clock()>timeOutWarning) &&
                     (!triggeredTimeoutWarning)
                  ) {
                     tarch::parallel::Node::getInstance().writeTimeOutWarning(
                     "exahype::records::State",
                     "send(int)", destination,tag,1
                     );
                     triggeredTimeoutWarning = true;
                  }
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                     (clock()>timeOutShutdown)
                  ) {
                     tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                     "exahype::records::State",
                     "send(int)", destination,tag,1
                     );
                  }
                  
               tarch::parallel::Node::getInstance().receiveDanglingMessages();
               usleep(communicateSleep);
               }
               
               delete sendRequestHandle;
               #ifdef Debug
               _log.debug("send(int,int)", "sent " + toString() );
               #endif
               
            }
            
         }
         
         
         
         void exahype::records::State::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
            if (communicateSleep<0) {
            
               MPI_Status  status;
               const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
               _senderDestinationRank = status.MPI_SOURCE;
               if ( result != MPI_SUCCESS ) {
                  std::ostringstream msg;
                  msg << "failed to start to receive exahype::records::State from node "
                  << source << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "receive(int)", msg.str() );
               }
               
            }
            else {
            
               MPI_Request* sendRequestHandle = new MPI_Request();
               MPI_Status   status;
               int          flag = 0;
               int          result;
               
               clock_t      timeOutWarning   = -1;
               clock_t      timeOutShutdown  = -1;
               bool         triggeredTimeoutWarning = false;
               
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  result = MPI_Irecv(
                     this, 1, Datatype, source, tag,
                     tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
                  );
                  
               }
               else {
                  result = MPI_Irecv(
                     this, 1, FullDatatype, source, tag,
                     tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
                  );
                  
               }
               if ( result != MPI_SUCCESS ) {
                  std::ostringstream msg;
                  msg << "failed to start to receive exahype::records::State from node "
                  << source << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "receive(int)", msg.str() );
               }
               
               result = MPI_Test( sendRequestHandle, &flag, &status );
               while (!flag) {
                  if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
                  if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
                  result = MPI_Test( sendRequestHandle, &flag, &status );
                  if (result!=MPI_SUCCESS) {
                     std::ostringstream msg;
                     msg << "testing for finished receive task for exahype::records::State failed: "
                     << tarch::parallel::MPIReturnValueToString(result);
                     _log.error("receive(int)", msg.str() );
                  }
                  
                  // deadlock aspect
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                     (clock()>timeOutWarning) &&
                     (!triggeredTimeoutWarning)
                  ) {
                     tarch::parallel::Node::getInstance().writeTimeOutWarning(
                     "exahype::records::State",
                     "receive(int)", source,tag,1
                     );
                     triggeredTimeoutWarning = true;
                  }
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                     (clock()>timeOutShutdown)
                  ) {
                     tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                     "exahype::records::State",
                     "receive(int)", source,tag,1
                     );
                  }
                  tarch::parallel::Node::getInstance().receiveDanglingMessages();
                  usleep(communicateSleep);
                  
               }
               
               delete sendRequestHandle;
               
               _senderDestinationRank = status.MPI_SOURCE;
               #ifdef Debug
               _log.debug("receive(int,int)", "received " + toString() ); 
               #endif
               
            }
            
         }
         
         
         
         bool exahype::records::State::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
            MPI_Status status;
            int  flag        = 0;
            MPI_Iprobe(
               MPI_ANY_SOURCE, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
            );
            if (flag) {
               int  messageCounter;
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  MPI_Get_count(&status, Datatype, &messageCounter);
               }
               else {
                  MPI_Get_count(&status, FullDatatype, &messageCounter);
               }
               return messageCounter > 0;
            }
            else return false;
            
         }
         
         int exahype::records::State::getSenderRank() const {
            assertion( _senderDestinationRank!=-1 );
            return _senderDestinationRank;
            
         }
      #endif
      
      
      exahype::records::StatePacked::PersistentRecords::PersistentRecords() {
         if ((9 >= (8 * sizeof(short int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((9 < (8 * sizeof(short int))));
         
      }
      
      
      exahype::records::StatePacked::PersistentRecords::PersistentRecords(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
      _maxRefinementLevelAllowed(maxRefinementLevelAllowed),
      _mergeMode(mergeMode),
      _sendMode(sendMode),
      _algorithmSection(algorithmSection),
      _minMeshWidth(minMeshWidth),
      _maxMeshWidth(maxMeshWidth),
      _numberOfInnerVertices(numberOfInnerVertices),
      _numberOfBoundaryVertices(numberOfBoundaryVertices),
      _numberOfOuterVertices(numberOfOuterVertices),
      _numberOfInnerCells(numberOfInnerCells),
      _numberOfOuterCells(numberOfOuterCells),
      _numberOfInnerLeafVertices(numberOfInnerLeafVertices),
      _numberOfBoundaryLeafVertices(numberOfBoundaryLeafVertices),
      _numberOfOuterLeafVertices(numberOfOuterLeafVertices),
      _numberOfInnerLeafCells(numberOfInnerLeafCells),
      _numberOfOuterLeafCells(numberOfOuterLeafCells),
      _maxLevel(maxLevel),
      _isTraversalInverted(isTraversalInverted) {
         setHasRefined(hasRefined);
         setHasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration);
         setHasErased(hasErased);
         setHasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration);
         setHasChangedVertexOrCellState(hasChangedVertexOrCellState);
         setHasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration);
         setReduceStateAndCell(reduceStateAndCell);
         setCouldNotEraseDueToDecompositionFlag(couldNotEraseDueToDecompositionFlag);
         setSubWorkerIsInvolvedInJoinOrFork(subWorkerIsInvolvedInJoinOrFork);
         if ((9 >= (8 * sizeof(short int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((9 < (8 * sizeof(short int))));
         
      }
      
      exahype::records::StatePacked::StatePacked() {
         if ((9 >= (8 * sizeof(short int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((9 < (8 * sizeof(short int))));
         
      }
      
      
      exahype::records::StatePacked::StatePacked(const PersistentRecords& persistentRecords):
      _persistentRecords(persistentRecords._maxRefinementLevelAllowed, persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._algorithmSection, persistentRecords._minMeshWidth, persistentRecords._maxMeshWidth, persistentRecords._numberOfInnerVertices, persistentRecords._numberOfBoundaryVertices, persistentRecords._numberOfOuterVertices, persistentRecords._numberOfInnerCells, persistentRecords._numberOfOuterCells, persistentRecords._numberOfInnerLeafVertices, persistentRecords._numberOfBoundaryLeafVertices, persistentRecords._numberOfOuterLeafVertices, persistentRecords._numberOfInnerLeafCells, persistentRecords._numberOfOuterLeafCells, persistentRecords._maxLevel, persistentRecords.getHasRefined(), persistentRecords.getHasTriggeredRefinementForNextIteration(), persistentRecords.getHasErased(), persistentRecords.getHasTriggeredEraseForNextIteration(), persistentRecords.getHasChangedVertexOrCellState(), persistentRecords.getHasModifiedGridInPreviousIteration(), persistentRecords._isTraversalInverted, persistentRecords.getReduceStateAndCell(), persistentRecords.getCouldNotEraseDueToDecompositionFlag(), persistentRecords.getSubWorkerIsInvolvedInJoinOrFork()) {
         if ((9 >= (8 * sizeof(short int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((9 < (8 * sizeof(short int))));
         
      }
      
      
      exahype::records::StatePacked::StatePacked(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
      _persistentRecords(maxRefinementLevelAllowed, mergeMode, sendMode, algorithmSection, minMeshWidth, maxMeshWidth, numberOfInnerVertices, numberOfBoundaryVertices, numberOfOuterVertices, numberOfInnerCells, numberOfOuterCells, numberOfInnerLeafVertices, numberOfBoundaryLeafVertices, numberOfOuterLeafVertices, numberOfInnerLeafCells, numberOfOuterLeafCells, maxLevel, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted, reduceStateAndCell, couldNotEraseDueToDecompositionFlag, subWorkerIsInvolvedInJoinOrFork) {
         if ((9 >= (8 * sizeof(short int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((9 < (8 * sizeof(short int))));
         
      }
      
      
      exahype::records::StatePacked::~StatePacked() { }
      
      std::string exahype::records::StatePacked::toString(const MergeMode& param) {
         return exahype::records::State::toString(param);
      }
      
      std::string exahype::records::StatePacked::getMergeModeMapping() {
         return exahype::records::State::getMergeModeMapping();
      }
      
      std::string exahype::records::StatePacked::toString(const SendMode& param) {
         return exahype::records::State::toString(param);
      }
      
      std::string exahype::records::StatePacked::getSendModeMapping() {
         return exahype::records::State::getSendModeMapping();
      }
      
      std::string exahype::records::StatePacked::toString(const AlgorithmSection& param) {
         return exahype::records::State::toString(param);
      }
      
      std::string exahype::records::StatePacked::getAlgorithmSectionMapping() {
         return exahype::records::State::getAlgorithmSectionMapping();
      }
      
      
      
      std::string exahype::records::StatePacked::toString() const {
         std::ostringstream stringstr;
         toString(stringstr);
         return stringstr.str();
      }
      
      void exahype::records::StatePacked::toString (std::ostream& out) const {
         out << "("; 
         out << "maxRefinementLevelAllowed:" << getMaxRefinementLevelAllowed();
         out << ",";
         out << "mergeMode:" << toString(getMergeMode());
         out << ",";
         out << "sendMode:" << toString(getSendMode());
         out << ",";
         out << "_algorithmSection:" << toString(getAlgorithmSection());
         out << ",";
         out << "minMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMinMeshWidth(i) << ",";
   }
   out << getMinMeshWidth(DIMENSIONS-1) << "]";
         out << ",";
         out << "maxMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMaxMeshWidth(i) << ",";
   }
   out << getMaxMeshWidth(DIMENSIONS-1) << "]";
         out << ",";
         out << "numberOfInnerVertices:" << getNumberOfInnerVertices();
         out << ",";
         out << "numberOfBoundaryVertices:" << getNumberOfBoundaryVertices();
         out << ",";
         out << "numberOfOuterVertices:" << getNumberOfOuterVertices();
         out << ",";
         out << "numberOfInnerCells:" << getNumberOfInnerCells();
         out << ",";
         out << "numberOfOuterCells:" << getNumberOfOuterCells();
         out << ",";
         out << "numberOfInnerLeafVertices:" << getNumberOfInnerLeafVertices();
         out << ",";
         out << "numberOfBoundaryLeafVertices:" << getNumberOfBoundaryLeafVertices();
         out << ",";
         out << "numberOfOuterLeafVertices:" << getNumberOfOuterLeafVertices();
         out << ",";
         out << "numberOfInnerLeafCells:" << getNumberOfInnerLeafCells();
         out << ",";
         out << "numberOfOuterLeafCells:" << getNumberOfOuterLeafCells();
         out << ",";
         out << "maxLevel:" << getMaxLevel();
         out << ",";
         out << "hasRefined:" << getHasRefined();
         out << ",";
         out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
         out << ",";
         out << "hasErased:" << getHasErased();
         out << ",";
         out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
         out << ",";
         out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
         out << ",";
         out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
         out << ",";
         out << "isTraversalInverted:" << getIsTraversalInverted();
         out << ",";
         out << "reduceStateAndCell:" << getReduceStateAndCell();
         out << ",";
         out << "couldNotEraseDueToDecompositionFlag:" << getCouldNotEraseDueToDecompositionFlag();
         out << ",";
         out << "subWorkerIsInvolvedInJoinOrFork:" << getSubWorkerIsInvolvedInJoinOrFork();
         out <<  ")";
      }
      
      
      exahype::records::StatePacked::PersistentRecords exahype::records::StatePacked::getPersistentRecords() const {
         return _persistentRecords;
      }
      
      exahype::records::State exahype::records::StatePacked::convert() const{
         return State(
            getMaxRefinementLevelAllowed(),
            getMergeMode(),
            getSendMode(),
            getAlgorithmSection(),
            getMinMeshWidth(),
            getMaxMeshWidth(),
            getNumberOfInnerVertices(),
            getNumberOfBoundaryVertices(),
            getNumberOfOuterVertices(),
            getNumberOfInnerCells(),
            getNumberOfOuterCells(),
            getNumberOfInnerLeafVertices(),
            getNumberOfBoundaryLeafVertices(),
            getNumberOfOuterLeafVertices(),
            getNumberOfInnerLeafCells(),
            getNumberOfOuterLeafCells(),
            getMaxLevel(),
            getHasRefined(),
            getHasTriggeredRefinementForNextIteration(),
            getHasErased(),
            getHasTriggeredEraseForNextIteration(),
            getHasChangedVertexOrCellState(),
            getHasModifiedGridInPreviousIteration(),
            getIsTraversalInverted(),
            getReduceStateAndCell(),
            getCouldNotEraseDueToDecompositionFlag(),
            getSubWorkerIsInvolvedInJoinOrFork()
         );
      }
      
      #ifdef Parallel
         tarch::logging::Log exahype::records::StatePacked::_log( "exahype::records::StatePacked" );
         
         MPI_Datatype exahype::records::StatePacked::Datatype = 0;
         MPI_Datatype exahype::records::StatePacked::FullDatatype = 0;
         
         
         void exahype::records::StatePacked::initDatatype() {
            {
               StatePacked dummyStatePacked[2];
               
               #ifdef MPI2
               const int Attributes = 19;
               #else
               const int Attributes = 20;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //maxRefinementLevelAllowed
                  , MPI_INT		 //mergeMode
                  , MPI_INT		 //sendMode
                  , MPI_INT		 //_algorithmSection
                  , MPI_DOUBLE		 //minMeshWidth
                  , MPI_DOUBLE		 //maxMeshWidth
                  , MPI_DOUBLE		 //numberOfInnerVertices
                  , MPI_DOUBLE		 //numberOfBoundaryVertices
                  , MPI_DOUBLE		 //numberOfOuterVertices
                  , MPI_DOUBLE		 //numberOfInnerCells
                  , MPI_DOUBLE		 //numberOfOuterCells
                  , MPI_DOUBLE		 //numberOfInnerLeafVertices
                  , MPI_DOUBLE		 //numberOfBoundaryLeafVertices
                  , MPI_DOUBLE		 //numberOfOuterLeafVertices
                  , MPI_DOUBLE		 //numberOfInnerLeafCells
                  , MPI_DOUBLE		 //numberOfOuterLeafCells
                  , MPI_INT		 //maxLevel
                  , MPI_CXX_BOOL		 //isTraversalInverted
                  , MPI_SHORT		 //_packedRecords0
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //maxRefinementLevelAllowed
                  , 1		 //mergeMode
                  , 1		 //sendMode
                  , 1		 //_algorithmSection
                  , DIMENSIONS		 //minMeshWidth
                  , DIMENSIONS		 //maxMeshWidth
                  , 1		 //numberOfInnerVertices
                  , 1		 //numberOfBoundaryVertices
                  , 1		 //numberOfOuterVertices
                  , 1		 //numberOfInnerCells
                  , 1		 //numberOfOuterCells
                  , 1		 //numberOfInnerLeafVertices
                  , 1		 //numberOfBoundaryLeafVertices
                  , 1		 //numberOfOuterLeafVertices
                  , 1		 //numberOfInnerLeafCells
                  , 1		 //numberOfOuterLeafCells
                  , 1		 //maxLevel
                  , 1		 //isTraversalInverted
                  , 1		 //_packedRecords0
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[18] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(StatePacked)), i, disp[i], Attributes, sizeof(StatePacked));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]))), 		&disp[19] );
               disp[19] -= base;
               disp[19] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &StatePacked::Datatype );
               MPI_Type_commit( &StatePacked::Datatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::Datatype);
               MPI_Type_commit( &StatePacked::Datatype );
               #endif
               
            }
            {
               StatePacked dummyStatePacked[2];
               
               #ifdef MPI2
               const int Attributes = 19;
               #else
               const int Attributes = 20;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //maxRefinementLevelAllowed
                  , MPI_INT		 //mergeMode
                  , MPI_INT		 //sendMode
                  , MPI_INT		 //_algorithmSection
                  , MPI_DOUBLE		 //minMeshWidth
                  , MPI_DOUBLE		 //maxMeshWidth
                  , MPI_DOUBLE		 //numberOfInnerVertices
                  , MPI_DOUBLE		 //numberOfBoundaryVertices
                  , MPI_DOUBLE		 //numberOfOuterVertices
                  , MPI_DOUBLE		 //numberOfInnerCells
                  , MPI_DOUBLE		 //numberOfOuterCells
                  , MPI_DOUBLE		 //numberOfInnerLeafVertices
                  , MPI_DOUBLE		 //numberOfBoundaryLeafVertices
                  , MPI_DOUBLE		 //numberOfOuterLeafVertices
                  , MPI_DOUBLE		 //numberOfInnerLeafCells
                  , MPI_DOUBLE		 //numberOfOuterLeafCells
                  , MPI_INT		 //maxLevel
                  , MPI_CXX_BOOL		 //isTraversalInverted
                  , MPI_SHORT		 //_packedRecords0
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //maxRefinementLevelAllowed
                  , 1		 //mergeMode
                  , 1		 //sendMode
                  , 1		 //_algorithmSection
                  , DIMENSIONS		 //minMeshWidth
                  , DIMENSIONS		 //maxMeshWidth
                  , 1		 //numberOfInnerVertices
                  , 1		 //numberOfBoundaryVertices
                  , 1		 //numberOfOuterVertices
                  , 1		 //numberOfInnerCells
                  , 1		 //numberOfOuterCells
                  , 1		 //numberOfInnerLeafVertices
                  , 1		 //numberOfBoundaryLeafVertices
                  , 1		 //numberOfOuterLeafVertices
                  , 1		 //numberOfInnerLeafCells
                  , 1		 //numberOfOuterLeafCells
                  , 1		 //maxLevel
                  , 1		 //isTraversalInverted
                  , 1		 //_packedRecords0
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[18] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(StatePacked)), i, disp[i], Attributes, sizeof(StatePacked));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]))), 		&disp[19] );
               disp[19] -= base;
               disp[19] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &StatePacked::FullDatatype );
               MPI_Type_commit( &StatePacked::FullDatatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::FullDatatype);
               MPI_Type_commit( &StatePacked::FullDatatype );
               #endif
               
            }
            
         }
         
         
         void exahype::records::StatePacked::shutdownDatatype() {
            MPI_Type_free( &StatePacked::Datatype );
            MPI_Type_free( &StatePacked::FullDatatype );
            
         }
         
         void exahype::records::StatePacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
            _senderDestinationRank = destination;
            
            if (communicateSleep<0) {
            
               const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
               if  (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "was not able to send message exahype::records::StatePacked "
                  << toString()
                  << " to node " << destination
                  << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "send(int)",msg.str() );
               }
               
            }
            else {
            
               MPI_Request* sendRequestHandle = new MPI_Request();
               MPI_Status   status;
               int          flag = 0;
               int          result;
               
               clock_t      timeOutWarning   = -1;
               clock_t      timeOutShutdown  = -1;
               bool         triggeredTimeoutWarning = false;
               
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  result = MPI_Isend(
                     this, 1, Datatype, destination,
                     tag, tarch::parallel::Node::getInstance().getCommunicator(),
                     sendRequestHandle
                  );
                  
               }
               else {
                  result = MPI_Isend(
                     this, 1, FullDatatype, destination,
                     tag, tarch::parallel::Node::getInstance().getCommunicator(),
                     sendRequestHandle
                  );
                  
               }
               if  (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "was not able to send message exahype::records::StatePacked "
                  << toString()
                  << " to node " << destination
                  << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "send(int)",msg.str() );
               }
               result = MPI_Test( sendRequestHandle, &flag, &status );
               while (!flag) {
                  if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
                  if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
                  result = MPI_Test( sendRequestHandle, &flag, &status );
                  if (result!=MPI_SUCCESS) {
                     std::ostringstream msg;
                     msg << "testing for finished send task for exahype::records::StatePacked "
                     << toString()
                     << " sent to node " << destination
                     << " failed: " << tarch::parallel::MPIReturnValueToString(result);
                     _log.error("send(int)", msg.str() );
                  }
                  
                  // deadlock aspect
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                     (clock()>timeOutWarning) &&
                     (!triggeredTimeoutWarning)
                  ) {
                     tarch::parallel::Node::getInstance().writeTimeOutWarning(
                     "exahype::records::StatePacked",
                     "send(int)", destination,tag,1
                     );
                     triggeredTimeoutWarning = true;
                  }
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                     (clock()>timeOutShutdown)
                  ) {
                     tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                     "exahype::records::StatePacked",
                     "send(int)", destination,tag,1
                     );
                  }
                  
               tarch::parallel::Node::getInstance().receiveDanglingMessages();
               usleep(communicateSleep);
               }
               
               delete sendRequestHandle;
               #ifdef Debug
               _log.debug("send(int,int)", "sent " + toString() );
               #endif
               
            }
            
         }
         
         
         
         void exahype::records::StatePacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
            if (communicateSleep<0) {
            
               MPI_Status  status;
               const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
               _senderDestinationRank = status.MPI_SOURCE;
               if ( result != MPI_SUCCESS ) {
                  std::ostringstream msg;
                  msg << "failed to start to receive exahype::records::StatePacked from node "
                  << source << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "receive(int)", msg.str() );
               }
               
            }
            else {
            
               MPI_Request* sendRequestHandle = new MPI_Request();
               MPI_Status   status;
               int          flag = 0;
               int          result;
               
               clock_t      timeOutWarning   = -1;
               clock_t      timeOutShutdown  = -1;
               bool         triggeredTimeoutWarning = false;
               
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  result = MPI_Irecv(
                     this, 1, Datatype, source, tag,
                     tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
                  );
                  
               }
               else {
                  result = MPI_Irecv(
                     this, 1, FullDatatype, source, tag,
                     tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
                  );
                  
               }
               if ( result != MPI_SUCCESS ) {
                  std::ostringstream msg;
                  msg << "failed to start to receive exahype::records::StatePacked from node "
                  << source << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "receive(int)", msg.str() );
               }
               
               result = MPI_Test( sendRequestHandle, &flag, &status );
               while (!flag) {
                  if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
                  if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
                  result = MPI_Test( sendRequestHandle, &flag, &status );
                  if (result!=MPI_SUCCESS) {
                     std::ostringstream msg;
                     msg << "testing for finished receive task for exahype::records::StatePacked failed: "
                     << tarch::parallel::MPIReturnValueToString(result);
                     _log.error("receive(int)", msg.str() );
                  }
                  
                  // deadlock aspect
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                     (clock()>timeOutWarning) &&
                     (!triggeredTimeoutWarning)
                  ) {
                     tarch::parallel::Node::getInstance().writeTimeOutWarning(
                     "exahype::records::StatePacked",
                     "receive(int)", source,tag,1
                     );
                     triggeredTimeoutWarning = true;
                  }
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                     (clock()>timeOutShutdown)
                  ) {
                     tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                     "exahype::records::StatePacked",
                     "receive(int)", source,tag,1
                     );
                  }
                  tarch::parallel::Node::getInstance().receiveDanglingMessages();
                  usleep(communicateSleep);
                  
               }
               
               delete sendRequestHandle;
               
               _senderDestinationRank = status.MPI_SOURCE;
               #ifdef Debug
               _log.debug("receive(int,int)", "received " + toString() ); 
               #endif
               
            }
            
         }
         
         
         
         bool exahype::records::StatePacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
            MPI_Status status;
            int  flag        = 0;
            MPI_Iprobe(
               MPI_ANY_SOURCE, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
            );
            if (flag) {
               int  messageCounter;
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  MPI_Get_count(&status, Datatype, &messageCounter);
               }
               else {
                  MPI_Get_count(&status, FullDatatype, &messageCounter);
               }
               return messageCounter > 0;
            }
            else return false;
            
         }
         
         int exahype::records::StatePacked::getSenderRank() const {
            assertion( _senderDestinationRank!=-1 );
            return _senderDestinationRank;
            
         }
      #endif
      
      
      
   #elif defined(TrackGridStatistics) && !defined(Parallel)
      exahype::records::State::PersistentRecords::PersistentRecords() {
         
      }
      
      
      exahype::records::State::PersistentRecords::PersistentRecords(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
      _maxRefinementLevelAllowed(maxRefinementLevelAllowed),
      _mergeMode(mergeMode),
      _sendMode(sendMode),
      _algorithmSection(algorithmSection),
      _minMeshWidth(minMeshWidth),
      _maxMeshWidth(maxMeshWidth),
      _numberOfInnerVertices(numberOfInnerVertices),
      _numberOfBoundaryVertices(numberOfBoundaryVertices),
      _numberOfOuterVertices(numberOfOuterVertices),
      _numberOfInnerCells(numberOfInnerCells),
      _numberOfOuterCells(numberOfOuterCells),
      _numberOfInnerLeafVertices(numberOfInnerLeafVertices),
      _numberOfBoundaryLeafVertices(numberOfBoundaryLeafVertices),
      _numberOfOuterLeafVertices(numberOfOuterLeafVertices),
      _numberOfInnerLeafCells(numberOfInnerLeafCells),
      _numberOfOuterLeafCells(numberOfOuterLeafCells),
      _maxLevel(maxLevel),
      _hasRefined(hasRefined),
      _hasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration),
      _hasErased(hasErased),
      _hasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration),
      _hasChangedVertexOrCellState(hasChangedVertexOrCellState),
      _hasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration),
      _isTraversalInverted(isTraversalInverted) {
         
      }
      
      exahype::records::State::State() {
         
      }
      
      
      exahype::records::State::State(const PersistentRecords& persistentRecords):
      _persistentRecords(persistentRecords._maxRefinementLevelAllowed, persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._algorithmSection, persistentRecords._minMeshWidth, persistentRecords._maxMeshWidth, persistentRecords._numberOfInnerVertices, persistentRecords._numberOfBoundaryVertices, persistentRecords._numberOfOuterVertices, persistentRecords._numberOfInnerCells, persistentRecords._numberOfOuterCells, persistentRecords._numberOfInnerLeafVertices, persistentRecords._numberOfBoundaryLeafVertices, persistentRecords._numberOfOuterLeafVertices, persistentRecords._numberOfInnerLeafCells, persistentRecords._numberOfOuterLeafCells, persistentRecords._maxLevel, persistentRecords._hasRefined, persistentRecords._hasTriggeredRefinementForNextIteration, persistentRecords._hasErased, persistentRecords._hasTriggeredEraseForNextIteration, persistentRecords._hasChangedVertexOrCellState, persistentRecords._hasModifiedGridInPreviousIteration, persistentRecords._isTraversalInverted) {
         
      }
      
      
      exahype::records::State::State(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
      _persistentRecords(maxRefinementLevelAllowed, mergeMode, sendMode, algorithmSection, minMeshWidth, maxMeshWidth, numberOfInnerVertices, numberOfBoundaryVertices, numberOfOuterVertices, numberOfInnerCells, numberOfOuterCells, numberOfInnerLeafVertices, numberOfBoundaryLeafVertices, numberOfOuterLeafVertices, numberOfInnerLeafCells, numberOfOuterLeafCells, maxLevel, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted) {
         
      }
      
      
      exahype::records::State::~State() { }
      
      std::string exahype::records::State::toString(const AlgorithmSection& param) {
         switch (param) {
            case TimeStepping: return "TimeStepping";
            case LimiterStatusSpreading: return "LimiterStatusSpreading";
            case MeshRefinement: return "MeshRefinement";
            case MeshRefinementOrLocalOrGlobalRecomputation: return "MeshRefinementOrLocalOrGlobalRecomputation";
            case LocalRecomputationAllSend: return "LocalRecomputationAllSend";
            case MeshRefinementOrGlobalRecomputation: return "MeshRefinementOrGlobalRecomputation";
            case MeshRefinementOrGlobalRecomputationAllSend: return "MeshRefinementOrGlobalRecomputationAllSend";
            case PredictionRerunAllSend: return "PredictionRerunAllSend";
         }
         return "undefined";
      }
      
      std::string exahype::records::State::getAlgorithmSectionMapping() {
         return "AlgorithmSection(TimeStepping=0,LimiterStatusSpreading=1,MeshRefinement=2,MeshRefinementOrLocalOrGlobalRecomputation=3,LocalRecomputationAllSend=4,MeshRefinementOrGlobalRecomputation=5,MeshRefinementOrGlobalRecomputationAllSend=6,PredictionRerunAllSend=7)";
      }
      std::string exahype::records::State::toString(const MergeMode& param) {
         switch (param) {
            case MergeNothing: return "MergeNothing";
            case BroadcastAndMergeTimeStepData: return "BroadcastAndMergeTimeStepData";
            case MergeFaceData: return "MergeFaceData";
            case DropFaceData: return "DropFaceData";
            case BroadcastAndMergeTimeStepDataAndMergeFaceData: return "BroadcastAndMergeTimeStepDataAndMergeFaceData";
            case BroadcastAndMergeTimeStepDataAndDropFaceData: return "BroadcastAndMergeTimeStepDataAndDropFaceData";
         }
         return "undefined";
      }
      
      std::string exahype::records::State::getMergeModeMapping() {
         return "MergeMode(MergeNothing=0,BroadcastAndMergeTimeStepData=1,MergeFaceData=2,DropFaceData=3,BroadcastAndMergeTimeStepDataAndMergeFaceData=4,BroadcastAndMergeTimeStepDataAndDropFaceData=5)";
      }
      std::string exahype::records::State::toString(const SendMode& param) {
         switch (param) {
            case SendNothing: return "SendNothing";
            case ReduceAndMergeTimeStepData: return "ReduceAndMergeTimeStepData";
            case SendFaceData: return "SendFaceData";
            case ReduceAndMergeTimeStepDataAndSendFaceData: return "ReduceAndMergeTimeStepDataAndSendFaceData";
         }
         return "undefined";
      }
      
      std::string exahype::records::State::getSendModeMapping() {
         return "SendMode(SendNothing=0,ReduceAndMergeTimeStepData=1,SendFaceData=2,ReduceAndMergeTimeStepDataAndSendFaceData=3)";
      }
      
      
      std::string exahype::records::State::toString() const {
         std::ostringstream stringstr;
         toString(stringstr);
         return stringstr.str();
      }
      
      void exahype::records::State::toString (std::ostream& out) const {
         out << "("; 
         out << "maxRefinementLevelAllowed:" << getMaxRefinementLevelAllowed();
         out << ",";
         out << "mergeMode:" << toString(getMergeMode());
         out << ",";
         out << "sendMode:" << toString(getSendMode());
         out << ",";
         out << "_algorithmSection:" << toString(getAlgorithmSection());
         out << ",";
         out << "minMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMinMeshWidth(i) << ",";
   }
   out << getMinMeshWidth(DIMENSIONS-1) << "]";
         out << ",";
         out << "maxMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMaxMeshWidth(i) << ",";
   }
   out << getMaxMeshWidth(DIMENSIONS-1) << "]";
         out << ",";
         out << "numberOfInnerVertices:" << getNumberOfInnerVertices();
         out << ",";
         out << "numberOfBoundaryVertices:" << getNumberOfBoundaryVertices();
         out << ",";
         out << "numberOfOuterVertices:" << getNumberOfOuterVertices();
         out << ",";
         out << "numberOfInnerCells:" << getNumberOfInnerCells();
         out << ",";
         out << "numberOfOuterCells:" << getNumberOfOuterCells();
         out << ",";
         out << "numberOfInnerLeafVertices:" << getNumberOfInnerLeafVertices();
         out << ",";
         out << "numberOfBoundaryLeafVertices:" << getNumberOfBoundaryLeafVertices();
         out << ",";
         out << "numberOfOuterLeafVertices:" << getNumberOfOuterLeafVertices();
         out << ",";
         out << "numberOfInnerLeafCells:" << getNumberOfInnerLeafCells();
         out << ",";
         out << "numberOfOuterLeafCells:" << getNumberOfOuterLeafCells();
         out << ",";
         out << "maxLevel:" << getMaxLevel();
         out << ",";
         out << "hasRefined:" << getHasRefined();
         out << ",";
         out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
         out << ",";
         out << "hasErased:" << getHasErased();
         out << ",";
         out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
         out << ",";
         out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
         out << ",";
         out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
         out << ",";
         out << "isTraversalInverted:" << getIsTraversalInverted();
         out <<  ")";
      }
      
      
      exahype::records::State::PersistentRecords exahype::records::State::getPersistentRecords() const {
         return _persistentRecords;
      }
      
      exahype::records::StatePacked exahype::records::State::convert() const{
         return StatePacked(
            getMaxRefinementLevelAllowed(),
            getMergeMode(),
            getSendMode(),
            getAlgorithmSection(),
            getMinMeshWidth(),
            getMaxMeshWidth(),
            getNumberOfInnerVertices(),
            getNumberOfBoundaryVertices(),
            getNumberOfOuterVertices(),
            getNumberOfInnerCells(),
            getNumberOfOuterCells(),
            getNumberOfInnerLeafVertices(),
            getNumberOfBoundaryLeafVertices(),
            getNumberOfOuterLeafVertices(),
            getNumberOfInnerLeafCells(),
            getNumberOfOuterLeafCells(),
            getMaxLevel(),
            getHasRefined(),
            getHasTriggeredRefinementForNextIteration(),
            getHasErased(),
            getHasTriggeredEraseForNextIteration(),
            getHasChangedVertexOrCellState(),
            getHasModifiedGridInPreviousIteration(),
            getIsTraversalInverted()
         );
      }
      
      #ifdef Parallel
         tarch::logging::Log exahype::records::State::_log( "exahype::records::State" );
         
         MPI_Datatype exahype::records::State::Datatype = 0;
         MPI_Datatype exahype::records::State::FullDatatype = 0;
         
         
         void exahype::records::State::initDatatype() {
            {
               State dummyState[2];
               
               #ifdef MPI2
               const int Attributes = 24;
               #else
               const int Attributes = 25;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //maxRefinementLevelAllowed
                  , MPI_INT		 //mergeMode
                  , MPI_INT		 //sendMode
                  , MPI_INT		 //_algorithmSection
                  , MPI_DOUBLE		 //minMeshWidth
                  , MPI_DOUBLE		 //maxMeshWidth
                  , MPI_DOUBLE		 //numberOfInnerVertices
                  , MPI_DOUBLE		 //numberOfBoundaryVertices
                  , MPI_DOUBLE		 //numberOfOuterVertices
                  , MPI_DOUBLE		 //numberOfInnerCells
                  , MPI_DOUBLE		 //numberOfOuterCells
                  , MPI_DOUBLE		 //numberOfInnerLeafVertices
                  , MPI_DOUBLE		 //numberOfBoundaryLeafVertices
                  , MPI_DOUBLE		 //numberOfOuterLeafVertices
                  , MPI_DOUBLE		 //numberOfInnerLeafCells
                  , MPI_DOUBLE		 //numberOfOuterLeafCells
                  , MPI_INT		 //maxLevel
                  , MPI_CXX_BOOL		 //hasRefined
                  , MPI_CXX_BOOL		 //hasTriggeredRefinementForNextIteration
                  , MPI_CXX_BOOL		 //hasErased
                  , MPI_CXX_BOOL		 //hasTriggeredEraseForNextIteration
                  , MPI_CXX_BOOL		 //hasChangedVertexOrCellState
                  , MPI_CXX_BOOL		 //hasModifiedGridInPreviousIteration
                  , MPI_CXX_BOOL		 //isTraversalInverted
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //maxRefinementLevelAllowed
                  , 1		 //mergeMode
                  , 1		 //sendMode
                  , 1		 //_algorithmSection
                  , DIMENSIONS		 //minMeshWidth
                  , DIMENSIONS		 //maxMeshWidth
                  , 1		 //numberOfInnerVertices
                  , 1		 //numberOfBoundaryVertices
                  , 1		 //numberOfOuterVertices
                  , 1		 //numberOfInnerCells
                  , 1		 //numberOfOuterCells
                  , 1		 //numberOfInnerLeafVertices
                  , 1		 //numberOfBoundaryLeafVertices
                  , 1		 //numberOfOuterLeafVertices
                  , 1		 //numberOfInnerLeafCells
                  , 1		 //numberOfOuterLeafCells
                  , 1		 //maxLevel
                  , 1		 //hasRefined
                  , 1		 //hasTriggeredRefinementForNextIteration
                  , 1		 //hasErased
                  , 1		 //hasTriggeredEraseForNextIteration
                  , 1		 //hasChangedVertexOrCellState
                  , 1		 //hasModifiedGridInPreviousIteration
                  , 1		 //isTraversalInverted
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[18] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[19] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[19] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[20] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[20] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[21] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[21] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[22] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[22] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[23] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[23] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(State)), i, disp[i], Attributes, sizeof(State));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]))), 		&disp[24] );
               disp[24] -= base;
               disp[24] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &State::Datatype );
               MPI_Type_commit( &State::Datatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::Datatype);
               MPI_Type_commit( &State::Datatype );
               #endif
               
            }
            {
               State dummyState[2];
               
               #ifdef MPI2
               const int Attributes = 24;
               #else
               const int Attributes = 25;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //maxRefinementLevelAllowed
                  , MPI_INT		 //mergeMode
                  , MPI_INT		 //sendMode
                  , MPI_INT		 //_algorithmSection
                  , MPI_DOUBLE		 //minMeshWidth
                  , MPI_DOUBLE		 //maxMeshWidth
                  , MPI_DOUBLE		 //numberOfInnerVertices
                  , MPI_DOUBLE		 //numberOfBoundaryVertices
                  , MPI_DOUBLE		 //numberOfOuterVertices
                  , MPI_DOUBLE		 //numberOfInnerCells
                  , MPI_DOUBLE		 //numberOfOuterCells
                  , MPI_DOUBLE		 //numberOfInnerLeafVertices
                  , MPI_DOUBLE		 //numberOfBoundaryLeafVertices
                  , MPI_DOUBLE		 //numberOfOuterLeafVertices
                  , MPI_DOUBLE		 //numberOfInnerLeafCells
                  , MPI_DOUBLE		 //numberOfOuterLeafCells
                  , MPI_INT		 //maxLevel
                  , MPI_CXX_BOOL		 //hasRefined
                  , MPI_CXX_BOOL		 //hasTriggeredRefinementForNextIteration
                  , MPI_CXX_BOOL		 //hasErased
                  , MPI_CXX_BOOL		 //hasTriggeredEraseForNextIteration
                  , MPI_CXX_BOOL		 //hasChangedVertexOrCellState
                  , MPI_CXX_BOOL		 //hasModifiedGridInPreviousIteration
                  , MPI_CXX_BOOL		 //isTraversalInverted
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //maxRefinementLevelAllowed
                  , 1		 //mergeMode
                  , 1		 //sendMode
                  , 1		 //_algorithmSection
                  , DIMENSIONS		 //minMeshWidth
                  , DIMENSIONS		 //maxMeshWidth
                  , 1		 //numberOfInnerVertices
                  , 1		 //numberOfBoundaryVertices
                  , 1		 //numberOfOuterVertices
                  , 1		 //numberOfInnerCells
                  , 1		 //numberOfOuterCells
                  , 1		 //numberOfInnerLeafVertices
                  , 1		 //numberOfBoundaryLeafVertices
                  , 1		 //numberOfOuterLeafVertices
                  , 1		 //numberOfInnerLeafCells
                  , 1		 //numberOfOuterLeafCells
                  , 1		 //maxLevel
                  , 1		 //hasRefined
                  , 1		 //hasTriggeredRefinementForNextIteration
                  , 1		 //hasErased
                  , 1		 //hasTriggeredEraseForNextIteration
                  , 1		 //hasChangedVertexOrCellState
                  , 1		 //hasModifiedGridInPreviousIteration
                  , 1		 //isTraversalInverted
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[18] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[19] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[19] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[20] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[20] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[21] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[21] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[22] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[22] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[23] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[23] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(State)), i, disp[i], Attributes, sizeof(State));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]))), 		&disp[24] );
               disp[24] -= base;
               disp[24] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &State::FullDatatype );
               MPI_Type_commit( &State::FullDatatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::FullDatatype);
               MPI_Type_commit( &State::FullDatatype );
               #endif
               
            }
            
         }
         
         
         void exahype::records::State::shutdownDatatype() {
            MPI_Type_free( &State::Datatype );
            MPI_Type_free( &State::FullDatatype );
            
         }
         
         void exahype::records::State::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
            _senderDestinationRank = destination;
            
            if (communicateSleep<0) {
            
               const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
               if  (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "was not able to send message exahype::records::State "
                  << toString()
                  << " to node " << destination
                  << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "send(int)",msg.str() );
               }
               
            }
            else {
            
               MPI_Request* sendRequestHandle = new MPI_Request();
               MPI_Status   status;
               int          flag = 0;
               int          result;
               
               clock_t      timeOutWarning   = -1;
               clock_t      timeOutShutdown  = -1;
               bool         triggeredTimeoutWarning = false;
               
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  result = MPI_Isend(
                     this, 1, Datatype, destination,
                     tag, tarch::parallel::Node::getInstance().getCommunicator(),
                     sendRequestHandle
                  );
                  
               }
               else {
                  result = MPI_Isend(
                     this, 1, FullDatatype, destination,
                     tag, tarch::parallel::Node::getInstance().getCommunicator(),
                     sendRequestHandle
                  );
                  
               }
               if  (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "was not able to send message exahype::records::State "
                  << toString()
                  << " to node " << destination
                  << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "send(int)",msg.str() );
               }
               result = MPI_Test( sendRequestHandle, &flag, &status );
               while (!flag) {
                  if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
                  if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
                  result = MPI_Test( sendRequestHandle, &flag, &status );
                  if (result!=MPI_SUCCESS) {
                     std::ostringstream msg;
                     msg << "testing for finished send task for exahype::records::State "
                     << toString()
                     << " sent to node " << destination
                     << " failed: " << tarch::parallel::MPIReturnValueToString(result);
                     _log.error("send(int)", msg.str() );
                  }
                  
                  // deadlock aspect
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                     (clock()>timeOutWarning) &&
                     (!triggeredTimeoutWarning)
                  ) {
                     tarch::parallel::Node::getInstance().writeTimeOutWarning(
                     "exahype::records::State",
                     "send(int)", destination,tag,1
                     );
                     triggeredTimeoutWarning = true;
                  }
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                     (clock()>timeOutShutdown)
                  ) {
                     tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                     "exahype::records::State",
                     "send(int)", destination,tag,1
                     );
                  }
                  
               tarch::parallel::Node::getInstance().receiveDanglingMessages();
               usleep(communicateSleep);
               }
               
               delete sendRequestHandle;
               #ifdef Debug
               _log.debug("send(int,int)", "sent " + toString() );
               #endif
               
            }
            
         }
         
         
         
         void exahype::records::State::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
            if (communicateSleep<0) {
            
               MPI_Status  status;
               const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
               _senderDestinationRank = status.MPI_SOURCE;
               if ( result != MPI_SUCCESS ) {
                  std::ostringstream msg;
                  msg << "failed to start to receive exahype::records::State from node "
                  << source << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "receive(int)", msg.str() );
               }
               
            }
            else {
            
               MPI_Request* sendRequestHandle = new MPI_Request();
               MPI_Status   status;
               int          flag = 0;
               int          result;
               
               clock_t      timeOutWarning   = -1;
               clock_t      timeOutShutdown  = -1;
               bool         triggeredTimeoutWarning = false;
               
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  result = MPI_Irecv(
                     this, 1, Datatype, source, tag,
                     tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
                  );
                  
               }
               else {
                  result = MPI_Irecv(
                     this, 1, FullDatatype, source, tag,
                     tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
                  );
                  
               }
               if ( result != MPI_SUCCESS ) {
                  std::ostringstream msg;
                  msg << "failed to start to receive exahype::records::State from node "
                  << source << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "receive(int)", msg.str() );
               }
               
               result = MPI_Test( sendRequestHandle, &flag, &status );
               while (!flag) {
                  if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
                  if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
                  result = MPI_Test( sendRequestHandle, &flag, &status );
                  if (result!=MPI_SUCCESS) {
                     std::ostringstream msg;
                     msg << "testing for finished receive task for exahype::records::State failed: "
                     << tarch::parallel::MPIReturnValueToString(result);
                     _log.error("receive(int)", msg.str() );
                  }
                  
                  // deadlock aspect
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                     (clock()>timeOutWarning) &&
                     (!triggeredTimeoutWarning)
                  ) {
                     tarch::parallel::Node::getInstance().writeTimeOutWarning(
                     "exahype::records::State",
                     "receive(int)", source,tag,1
                     );
                     triggeredTimeoutWarning = true;
                  }
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                     (clock()>timeOutShutdown)
                  ) {
                     tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                     "exahype::records::State",
                     "receive(int)", source,tag,1
                     );
                  }
                  tarch::parallel::Node::getInstance().receiveDanglingMessages();
                  usleep(communicateSleep);
                  
               }
               
               delete sendRequestHandle;
               
               _senderDestinationRank = status.MPI_SOURCE;
               #ifdef Debug
               _log.debug("receive(int,int)", "received " + toString() ); 
               #endif
               
            }
            
         }
         
         
         
         bool exahype::records::State::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
            MPI_Status status;
            int  flag        = 0;
            MPI_Iprobe(
               MPI_ANY_SOURCE, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
            );
            if (flag) {
               int  messageCounter;
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  MPI_Get_count(&status, Datatype, &messageCounter);
               }
               else {
                  MPI_Get_count(&status, FullDatatype, &messageCounter);
               }
               return messageCounter > 0;
            }
            else return false;
            
         }
         
         int exahype::records::State::getSenderRank() const {
            assertion( _senderDestinationRank!=-1 );
            return _senderDestinationRank;
            
         }
      #endif
      
      
      exahype::records::StatePacked::PersistentRecords::PersistentRecords() {
         if ((6 >= (8 * sizeof(short int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((6 < (8 * sizeof(short int))));
         
      }
      
      
      exahype::records::StatePacked::PersistentRecords::PersistentRecords(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
      _maxRefinementLevelAllowed(maxRefinementLevelAllowed),
      _mergeMode(mergeMode),
      _sendMode(sendMode),
      _algorithmSection(algorithmSection),
      _minMeshWidth(minMeshWidth),
      _maxMeshWidth(maxMeshWidth),
      _numberOfInnerVertices(numberOfInnerVertices),
      _numberOfBoundaryVertices(numberOfBoundaryVertices),
      _numberOfOuterVertices(numberOfOuterVertices),
      _numberOfInnerCells(numberOfInnerCells),
      _numberOfOuterCells(numberOfOuterCells),
      _numberOfInnerLeafVertices(numberOfInnerLeafVertices),
      _numberOfBoundaryLeafVertices(numberOfBoundaryLeafVertices),
      _numberOfOuterLeafVertices(numberOfOuterLeafVertices),
      _numberOfInnerLeafCells(numberOfInnerLeafCells),
      _numberOfOuterLeafCells(numberOfOuterLeafCells),
      _maxLevel(maxLevel),
      _isTraversalInverted(isTraversalInverted) {
         setHasRefined(hasRefined);
         setHasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration);
         setHasErased(hasErased);
         setHasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration);
         setHasChangedVertexOrCellState(hasChangedVertexOrCellState);
         setHasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration);
         if ((6 >= (8 * sizeof(short int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((6 < (8 * sizeof(short int))));
         
      }
      
      exahype::records::StatePacked::StatePacked() {
         if ((6 >= (8 * sizeof(short int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((6 < (8 * sizeof(short int))));
         
      }
      
      
      exahype::records::StatePacked::StatePacked(const PersistentRecords& persistentRecords):
      _persistentRecords(persistentRecords._maxRefinementLevelAllowed, persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._algorithmSection, persistentRecords._minMeshWidth, persistentRecords._maxMeshWidth, persistentRecords._numberOfInnerVertices, persistentRecords._numberOfBoundaryVertices, persistentRecords._numberOfOuterVertices, persistentRecords._numberOfInnerCells, persistentRecords._numberOfOuterCells, persistentRecords._numberOfInnerLeafVertices, persistentRecords._numberOfBoundaryLeafVertices, persistentRecords._numberOfOuterLeafVertices, persistentRecords._numberOfInnerLeafCells, persistentRecords._numberOfOuterLeafCells, persistentRecords._maxLevel, persistentRecords.getHasRefined(), persistentRecords.getHasTriggeredRefinementForNextIteration(), persistentRecords.getHasErased(), persistentRecords.getHasTriggeredEraseForNextIteration(), persistentRecords.getHasChangedVertexOrCellState(), persistentRecords.getHasModifiedGridInPreviousIteration(), persistentRecords._isTraversalInverted) {
         if ((6 >= (8 * sizeof(short int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((6 < (8 * sizeof(short int))));
         
      }
      
      
      exahype::records::StatePacked::StatePacked(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const tarch::la::Vector<DIMENSIONS,double>& minMeshWidth, const tarch::la::Vector<DIMENSIONS,double>& maxMeshWidth, const double& numberOfInnerVertices, const double& numberOfBoundaryVertices, const double& numberOfOuterVertices, const double& numberOfInnerCells, const double& numberOfOuterCells, const double& numberOfInnerLeafVertices, const double& numberOfBoundaryLeafVertices, const double& numberOfOuterLeafVertices, const double& numberOfInnerLeafCells, const double& numberOfOuterLeafCells, const int& maxLevel, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted):
      _persistentRecords(maxRefinementLevelAllowed, mergeMode, sendMode, algorithmSection, minMeshWidth, maxMeshWidth, numberOfInnerVertices, numberOfBoundaryVertices, numberOfOuterVertices, numberOfInnerCells, numberOfOuterCells, numberOfInnerLeafVertices, numberOfBoundaryLeafVertices, numberOfOuterLeafVertices, numberOfInnerLeafCells, numberOfOuterLeafCells, maxLevel, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted) {
         if ((6 >= (8 * sizeof(short int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((6 < (8 * sizeof(short int))));
         
      }
      
      
      exahype::records::StatePacked::~StatePacked() { }
      
      std::string exahype::records::StatePacked::toString(const MergeMode& param) {
         return exahype::records::State::toString(param);
      }
      
      std::string exahype::records::StatePacked::getMergeModeMapping() {
         return exahype::records::State::getMergeModeMapping();
      }
      
      std::string exahype::records::StatePacked::toString(const SendMode& param) {
         return exahype::records::State::toString(param);
      }
      
      std::string exahype::records::StatePacked::getSendModeMapping() {
         return exahype::records::State::getSendModeMapping();
      }
      
      std::string exahype::records::StatePacked::toString(const AlgorithmSection& param) {
         return exahype::records::State::toString(param);
      }
      
      std::string exahype::records::StatePacked::getAlgorithmSectionMapping() {
         return exahype::records::State::getAlgorithmSectionMapping();
      }
      
      
      
      std::string exahype::records::StatePacked::toString() const {
         std::ostringstream stringstr;
         toString(stringstr);
         return stringstr.str();
      }
      
      void exahype::records::StatePacked::toString (std::ostream& out) const {
         out << "("; 
         out << "maxRefinementLevelAllowed:" << getMaxRefinementLevelAllowed();
         out << ",";
         out << "mergeMode:" << toString(getMergeMode());
         out << ",";
         out << "sendMode:" << toString(getSendMode());
         out << ",";
         out << "_algorithmSection:" << toString(getAlgorithmSection());
         out << ",";
         out << "minMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMinMeshWidth(i) << ",";
   }
   out << getMinMeshWidth(DIMENSIONS-1) << "]";
         out << ",";
         out << "maxMeshWidth:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getMaxMeshWidth(i) << ",";
   }
   out << getMaxMeshWidth(DIMENSIONS-1) << "]";
         out << ",";
         out << "numberOfInnerVertices:" << getNumberOfInnerVertices();
         out << ",";
         out << "numberOfBoundaryVertices:" << getNumberOfBoundaryVertices();
         out << ",";
         out << "numberOfOuterVertices:" << getNumberOfOuterVertices();
         out << ",";
         out << "numberOfInnerCells:" << getNumberOfInnerCells();
         out << ",";
         out << "numberOfOuterCells:" << getNumberOfOuterCells();
         out << ",";
         out << "numberOfInnerLeafVertices:" << getNumberOfInnerLeafVertices();
         out << ",";
         out << "numberOfBoundaryLeafVertices:" << getNumberOfBoundaryLeafVertices();
         out << ",";
         out << "numberOfOuterLeafVertices:" << getNumberOfOuterLeafVertices();
         out << ",";
         out << "numberOfInnerLeafCells:" << getNumberOfInnerLeafCells();
         out << ",";
         out << "numberOfOuterLeafCells:" << getNumberOfOuterLeafCells();
         out << ",";
         out << "maxLevel:" << getMaxLevel();
         out << ",";
         out << "hasRefined:" << getHasRefined();
         out << ",";
         out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
         out << ",";
         out << "hasErased:" << getHasErased();
         out << ",";
         out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
         out << ",";
         out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
         out << ",";
         out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
         out << ",";
         out << "isTraversalInverted:" << getIsTraversalInverted();
         out <<  ")";
      }
      
      
      exahype::records::StatePacked::PersistentRecords exahype::records::StatePacked::getPersistentRecords() const {
         return _persistentRecords;
      }
      
      exahype::records::State exahype::records::StatePacked::convert() const{
         return State(
            getMaxRefinementLevelAllowed(),
            getMergeMode(),
            getSendMode(),
            getAlgorithmSection(),
            getMinMeshWidth(),
            getMaxMeshWidth(),
            getNumberOfInnerVertices(),
            getNumberOfBoundaryVertices(),
            getNumberOfOuterVertices(),
            getNumberOfInnerCells(),
            getNumberOfOuterCells(),
            getNumberOfInnerLeafVertices(),
            getNumberOfBoundaryLeafVertices(),
            getNumberOfOuterLeafVertices(),
            getNumberOfInnerLeafCells(),
            getNumberOfOuterLeafCells(),
            getMaxLevel(),
            getHasRefined(),
            getHasTriggeredRefinementForNextIteration(),
            getHasErased(),
            getHasTriggeredEraseForNextIteration(),
            getHasChangedVertexOrCellState(),
            getHasModifiedGridInPreviousIteration(),
            getIsTraversalInverted()
         );
      }
      
      #ifdef Parallel
         tarch::logging::Log exahype::records::StatePacked::_log( "exahype::records::StatePacked" );
         
         MPI_Datatype exahype::records::StatePacked::Datatype = 0;
         MPI_Datatype exahype::records::StatePacked::FullDatatype = 0;
         
         
         void exahype::records::StatePacked::initDatatype() {
            {
               StatePacked dummyStatePacked[2];
               
               #ifdef MPI2
               const int Attributes = 19;
               #else
               const int Attributes = 20;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //maxRefinementLevelAllowed
                  , MPI_INT		 //mergeMode
                  , MPI_INT		 //sendMode
                  , MPI_INT		 //_algorithmSection
                  , MPI_DOUBLE		 //minMeshWidth
                  , MPI_DOUBLE		 //maxMeshWidth
                  , MPI_DOUBLE		 //numberOfInnerVertices
                  , MPI_DOUBLE		 //numberOfBoundaryVertices
                  , MPI_DOUBLE		 //numberOfOuterVertices
                  , MPI_DOUBLE		 //numberOfInnerCells
                  , MPI_DOUBLE		 //numberOfOuterCells
                  , MPI_DOUBLE		 //numberOfInnerLeafVertices
                  , MPI_DOUBLE		 //numberOfBoundaryLeafVertices
                  , MPI_DOUBLE		 //numberOfOuterLeafVertices
                  , MPI_DOUBLE		 //numberOfInnerLeafCells
                  , MPI_DOUBLE		 //numberOfOuterLeafCells
                  , MPI_INT		 //maxLevel
                  , MPI_CXX_BOOL		 //isTraversalInverted
                  , MPI_SHORT		 //_packedRecords0
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //maxRefinementLevelAllowed
                  , 1		 //mergeMode
                  , 1		 //sendMode
                  , 1		 //_algorithmSection
                  , DIMENSIONS		 //minMeshWidth
                  , DIMENSIONS		 //maxMeshWidth
                  , 1		 //numberOfInnerVertices
                  , 1		 //numberOfBoundaryVertices
                  , 1		 //numberOfOuterVertices
                  , 1		 //numberOfInnerCells
                  , 1		 //numberOfOuterCells
                  , 1		 //numberOfInnerLeafVertices
                  , 1		 //numberOfBoundaryLeafVertices
                  , 1		 //numberOfOuterLeafVertices
                  , 1		 //numberOfInnerLeafCells
                  , 1		 //numberOfOuterLeafCells
                  , 1		 //maxLevel
                  , 1		 //isTraversalInverted
                  , 1		 //_packedRecords0
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[18] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(StatePacked)), i, disp[i], Attributes, sizeof(StatePacked));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]))), 		&disp[19] );
               disp[19] -= base;
               disp[19] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &StatePacked::Datatype );
               MPI_Type_commit( &StatePacked::Datatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::Datatype);
               MPI_Type_commit( &StatePacked::Datatype );
               #endif
               
            }
            {
               StatePacked dummyStatePacked[2];
               
               #ifdef MPI2
               const int Attributes = 19;
               #else
               const int Attributes = 20;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //maxRefinementLevelAllowed
                  , MPI_INT		 //mergeMode
                  , MPI_INT		 //sendMode
                  , MPI_INT		 //_algorithmSection
                  , MPI_DOUBLE		 //minMeshWidth
                  , MPI_DOUBLE		 //maxMeshWidth
                  , MPI_DOUBLE		 //numberOfInnerVertices
                  , MPI_DOUBLE		 //numberOfBoundaryVertices
                  , MPI_DOUBLE		 //numberOfOuterVertices
                  , MPI_DOUBLE		 //numberOfInnerCells
                  , MPI_DOUBLE		 //numberOfOuterCells
                  , MPI_DOUBLE		 //numberOfInnerLeafVertices
                  , MPI_DOUBLE		 //numberOfBoundaryLeafVertices
                  , MPI_DOUBLE		 //numberOfOuterLeafVertices
                  , MPI_DOUBLE		 //numberOfInnerLeafCells
                  , MPI_DOUBLE		 //numberOfOuterLeafCells
                  , MPI_INT		 //maxLevel
                  , MPI_CXX_BOOL		 //isTraversalInverted
                  , MPI_SHORT		 //_packedRecords0
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //maxRefinementLevelAllowed
                  , 1		 //mergeMode
                  , 1		 //sendMode
                  , 1		 //_algorithmSection
                  , DIMENSIONS		 //minMeshWidth
                  , DIMENSIONS		 //maxMeshWidth
                  , 1		 //numberOfInnerVertices
                  , 1		 //numberOfBoundaryVertices
                  , 1		 //numberOfOuterVertices
                  , 1		 //numberOfInnerCells
                  , 1		 //numberOfOuterCells
                  , 1		 //numberOfInnerLeafVertices
                  , 1		 //numberOfBoundaryLeafVertices
                  , 1		 //numberOfOuterLeafVertices
                  , 1		 //numberOfInnerLeafCells
                  , 1		 //numberOfOuterLeafCells
                  , 1		 //maxLevel
                  , 1		 //isTraversalInverted
                  , 1		 //_packedRecords0
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._minMeshWidth[0]))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxMeshWidth[0]))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerVertices))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryVertices))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterVertices))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerCells))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterCells))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafVertices))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfBoundaryLeafVertices))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafVertices))), 		&disp[13] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfInnerLeafCells))), 		&disp[14] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._numberOfOuterLeafCells))), 		&disp[15] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxLevel))), 		&disp[16] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[17] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[17] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[18] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[18] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(StatePacked)), i, disp[i], Attributes, sizeof(StatePacked));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]))), 		&disp[19] );
               disp[19] -= base;
               disp[19] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &StatePacked::FullDatatype );
               MPI_Type_commit( &StatePacked::FullDatatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::FullDatatype);
               MPI_Type_commit( &StatePacked::FullDatatype );
               #endif
               
            }
            
         }
         
         
         void exahype::records::StatePacked::shutdownDatatype() {
            MPI_Type_free( &StatePacked::Datatype );
            MPI_Type_free( &StatePacked::FullDatatype );
            
         }
         
         void exahype::records::StatePacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
            _senderDestinationRank = destination;
            
            if (communicateSleep<0) {
            
               const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
               if  (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "was not able to send message exahype::records::StatePacked "
                  << toString()
                  << " to node " << destination
                  << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "send(int)",msg.str() );
               }
               
            }
            else {
            
               MPI_Request* sendRequestHandle = new MPI_Request();
               MPI_Status   status;
               int          flag = 0;
               int          result;
               
               clock_t      timeOutWarning   = -1;
               clock_t      timeOutShutdown  = -1;
               bool         triggeredTimeoutWarning = false;
               
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  result = MPI_Isend(
                     this, 1, Datatype, destination,
                     tag, tarch::parallel::Node::getInstance().getCommunicator(),
                     sendRequestHandle
                  );
                  
               }
               else {
                  result = MPI_Isend(
                     this, 1, FullDatatype, destination,
                     tag, tarch::parallel::Node::getInstance().getCommunicator(),
                     sendRequestHandle
                  );
                  
               }
               if  (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "was not able to send message exahype::records::StatePacked "
                  << toString()
                  << " to node " << destination
                  << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "send(int)",msg.str() );
               }
               result = MPI_Test( sendRequestHandle, &flag, &status );
               while (!flag) {
                  if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
                  if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
                  result = MPI_Test( sendRequestHandle, &flag, &status );
                  if (result!=MPI_SUCCESS) {
                     std::ostringstream msg;
                     msg << "testing for finished send task for exahype::records::StatePacked "
                     << toString()
                     << " sent to node " << destination
                     << " failed: " << tarch::parallel::MPIReturnValueToString(result);
                     _log.error("send(int)", msg.str() );
                  }
                  
                  // deadlock aspect
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                     (clock()>timeOutWarning) &&
                     (!triggeredTimeoutWarning)
                  ) {
                     tarch::parallel::Node::getInstance().writeTimeOutWarning(
                     "exahype::records::StatePacked",
                     "send(int)", destination,tag,1
                     );
                     triggeredTimeoutWarning = true;
                  }
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                     (clock()>timeOutShutdown)
                  ) {
                     tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                     "exahype::records::StatePacked",
                     "send(int)", destination,tag,1
                     );
                  }
                  
               tarch::parallel::Node::getInstance().receiveDanglingMessages();
               usleep(communicateSleep);
               }
               
               delete sendRequestHandle;
               #ifdef Debug
               _log.debug("send(int,int)", "sent " + toString() );
               #endif
               
            }
            
         }
         
         
         
         void exahype::records::StatePacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
            if (communicateSleep<0) {
            
               MPI_Status  status;
               const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
               _senderDestinationRank = status.MPI_SOURCE;
               if ( result != MPI_SUCCESS ) {
                  std::ostringstream msg;
                  msg << "failed to start to receive exahype::records::StatePacked from node "
                  << source << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "receive(int)", msg.str() );
               }
               
            }
            else {
            
               MPI_Request* sendRequestHandle = new MPI_Request();
               MPI_Status   status;
               int          flag = 0;
               int          result;
               
               clock_t      timeOutWarning   = -1;
               clock_t      timeOutShutdown  = -1;
               bool         triggeredTimeoutWarning = false;
               
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  result = MPI_Irecv(
                     this, 1, Datatype, source, tag,
                     tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
                  );
                  
               }
               else {
                  result = MPI_Irecv(
                     this, 1, FullDatatype, source, tag,
                     tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
                  );
                  
               }
               if ( result != MPI_SUCCESS ) {
                  std::ostringstream msg;
                  msg << "failed to start to receive exahype::records::StatePacked from node "
                  << source << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "receive(int)", msg.str() );
               }
               
               result = MPI_Test( sendRequestHandle, &flag, &status );
               while (!flag) {
                  if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
                  if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
                  result = MPI_Test( sendRequestHandle, &flag, &status );
                  if (result!=MPI_SUCCESS) {
                     std::ostringstream msg;
                     msg << "testing for finished receive task for exahype::records::StatePacked failed: "
                     << tarch::parallel::MPIReturnValueToString(result);
                     _log.error("receive(int)", msg.str() );
                  }
                  
                  // deadlock aspect
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                     (clock()>timeOutWarning) &&
                     (!triggeredTimeoutWarning)
                  ) {
                     tarch::parallel::Node::getInstance().writeTimeOutWarning(
                     "exahype::records::StatePacked",
                     "receive(int)", source,tag,1
                     );
                     triggeredTimeoutWarning = true;
                  }
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                     (clock()>timeOutShutdown)
                  ) {
                     tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                     "exahype::records::StatePacked",
                     "receive(int)", source,tag,1
                     );
                  }
                  tarch::parallel::Node::getInstance().receiveDanglingMessages();
                  usleep(communicateSleep);
                  
               }
               
               delete sendRequestHandle;
               
               _senderDestinationRank = status.MPI_SOURCE;
               #ifdef Debug
               _log.debug("receive(int,int)", "received " + toString() ); 
               #endif
               
            }
            
         }
         
         
         
         bool exahype::records::StatePacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
            MPI_Status status;
            int  flag        = 0;
            MPI_Iprobe(
               MPI_ANY_SOURCE, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
            );
            if (flag) {
               int  messageCounter;
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  MPI_Get_count(&status, Datatype, &messageCounter);
               }
               else {
                  MPI_Get_count(&status, FullDatatype, &messageCounter);
               }
               return messageCounter > 0;
            }
            else return false;
            
         }
         
         int exahype::records::StatePacked::getSenderRank() const {
            assertion( _senderDestinationRank!=-1 );
            return _senderDestinationRank;
            
         }
      #endif
      
      
      
   #elif !defined(TrackGridStatistics) && defined(Parallel)
      exahype::records::State::PersistentRecords::PersistentRecords() {
         
      }
      
      
      exahype::records::State::PersistentRecords::PersistentRecords(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
      _maxRefinementLevelAllowed(maxRefinementLevelAllowed),
      _mergeMode(mergeMode),
      _sendMode(sendMode),
      _algorithmSection(algorithmSection),
      _hasRefined(hasRefined),
      _hasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration),
      _hasErased(hasErased),
      _hasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration),
      _hasChangedVertexOrCellState(hasChangedVertexOrCellState),
      _hasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration),
      _isTraversalInverted(isTraversalInverted),
      _reduceStateAndCell(reduceStateAndCell),
      _couldNotEraseDueToDecompositionFlag(couldNotEraseDueToDecompositionFlag),
      _subWorkerIsInvolvedInJoinOrFork(subWorkerIsInvolvedInJoinOrFork) {
         
      }
      
      exahype::records::State::State() {
         
      }
      
      
      exahype::records::State::State(const PersistentRecords& persistentRecords):
      _persistentRecords(persistentRecords._maxRefinementLevelAllowed, persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._algorithmSection, persistentRecords._hasRefined, persistentRecords._hasTriggeredRefinementForNextIteration, persistentRecords._hasErased, persistentRecords._hasTriggeredEraseForNextIteration, persistentRecords._hasChangedVertexOrCellState, persistentRecords._hasModifiedGridInPreviousIteration, persistentRecords._isTraversalInverted, persistentRecords._reduceStateAndCell, persistentRecords._couldNotEraseDueToDecompositionFlag, persistentRecords._subWorkerIsInvolvedInJoinOrFork) {
         
      }
      
      
      exahype::records::State::State(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
      _persistentRecords(maxRefinementLevelAllowed, mergeMode, sendMode, algorithmSection, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted, reduceStateAndCell, couldNotEraseDueToDecompositionFlag, subWorkerIsInvolvedInJoinOrFork) {
         
      }
      
      
      exahype::records::State::~State() { }
      
      std::string exahype::records::State::toString(const AlgorithmSection& param) {
         switch (param) {
            case TimeStepping: return "TimeStepping";
            case LimiterStatusSpreading: return "LimiterStatusSpreading";
            case MeshRefinement: return "MeshRefinement";
            case MeshRefinementOrLocalOrGlobalRecomputation: return "MeshRefinementOrLocalOrGlobalRecomputation";
            case LocalRecomputationAllSend: return "LocalRecomputationAllSend";
            case MeshRefinementOrGlobalRecomputation: return "MeshRefinementOrGlobalRecomputation";
            case MeshRefinementOrGlobalRecomputationAllSend: return "MeshRefinementOrGlobalRecomputationAllSend";
            case PredictionRerunAllSend: return "PredictionRerunAllSend";
         }
         return "undefined";
      }
      
      std::string exahype::records::State::getAlgorithmSectionMapping() {
         return "AlgorithmSection(TimeStepping=0,LimiterStatusSpreading=1,MeshRefinement=2,MeshRefinementOrLocalOrGlobalRecomputation=3,LocalRecomputationAllSend=4,MeshRefinementOrGlobalRecomputation=5,MeshRefinementOrGlobalRecomputationAllSend=6,PredictionRerunAllSend=7)";
      }
      std::string exahype::records::State::toString(const MergeMode& param) {
         switch (param) {
            case MergeNothing: return "MergeNothing";
            case BroadcastAndMergeTimeStepData: return "BroadcastAndMergeTimeStepData";
            case MergeFaceData: return "MergeFaceData";
            case DropFaceData: return "DropFaceData";
            case BroadcastAndMergeTimeStepDataAndMergeFaceData: return "BroadcastAndMergeTimeStepDataAndMergeFaceData";
            case BroadcastAndMergeTimeStepDataAndDropFaceData: return "BroadcastAndMergeTimeStepDataAndDropFaceData";
         }
         return "undefined";
      }
      
      std::string exahype::records::State::getMergeModeMapping() {
         return "MergeMode(MergeNothing=0,BroadcastAndMergeTimeStepData=1,MergeFaceData=2,DropFaceData=3,BroadcastAndMergeTimeStepDataAndMergeFaceData=4,BroadcastAndMergeTimeStepDataAndDropFaceData=5)";
      }
      std::string exahype::records::State::toString(const SendMode& param) {
         switch (param) {
            case SendNothing: return "SendNothing";
            case ReduceAndMergeTimeStepData: return "ReduceAndMergeTimeStepData";
            case SendFaceData: return "SendFaceData";
            case ReduceAndMergeTimeStepDataAndSendFaceData: return "ReduceAndMergeTimeStepDataAndSendFaceData";
         }
         return "undefined";
      }
      
      std::string exahype::records::State::getSendModeMapping() {
         return "SendMode(SendNothing=0,ReduceAndMergeTimeStepData=1,SendFaceData=2,ReduceAndMergeTimeStepDataAndSendFaceData=3)";
      }
      
      
      std::string exahype::records::State::toString() const {
         std::ostringstream stringstr;
         toString(stringstr);
         return stringstr.str();
      }
      
      void exahype::records::State::toString (std::ostream& out) const {
         out << "("; 
         out << "maxRefinementLevelAllowed:" << getMaxRefinementLevelAllowed();
         out << ",";
         out << "mergeMode:" << toString(getMergeMode());
         out << ",";
         out << "sendMode:" << toString(getSendMode());
         out << ",";
         out << "_algorithmSection:" << toString(getAlgorithmSection());
         out << ",";
         out << "hasRefined:" << getHasRefined();
         out << ",";
         out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
         out << ",";
         out << "hasErased:" << getHasErased();
         out << ",";
         out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
         out << ",";
         out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
         out << ",";
         out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
         out << ",";
         out << "isTraversalInverted:" << getIsTraversalInverted();
         out << ",";
         out << "reduceStateAndCell:" << getReduceStateAndCell();
         out << ",";
         out << "couldNotEraseDueToDecompositionFlag:" << getCouldNotEraseDueToDecompositionFlag();
         out << ",";
         out << "subWorkerIsInvolvedInJoinOrFork:" << getSubWorkerIsInvolvedInJoinOrFork();
         out <<  ")";
      }
      
      
      exahype::records::State::PersistentRecords exahype::records::State::getPersistentRecords() const {
         return _persistentRecords;
      }
      
      exahype::records::StatePacked exahype::records::State::convert() const{
         return StatePacked(
            getMaxRefinementLevelAllowed(),
            getMergeMode(),
            getSendMode(),
            getAlgorithmSection(),
            getHasRefined(),
            getHasTriggeredRefinementForNextIteration(),
            getHasErased(),
            getHasTriggeredEraseForNextIteration(),
            getHasChangedVertexOrCellState(),
            getHasModifiedGridInPreviousIteration(),
            getIsTraversalInverted(),
            getReduceStateAndCell(),
            getCouldNotEraseDueToDecompositionFlag(),
            getSubWorkerIsInvolvedInJoinOrFork()
         );
      }
      
      #ifdef Parallel
         tarch::logging::Log exahype::records::State::_log( "exahype::records::State" );
         
         MPI_Datatype exahype::records::State::Datatype = 0;
         MPI_Datatype exahype::records::State::FullDatatype = 0;
         
         
         void exahype::records::State::initDatatype() {
            {
               State dummyState[2];
               
               #ifdef MPI2
               const int Attributes = 14;
               #else
               const int Attributes = 15;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //maxRefinementLevelAllowed
                  , MPI_INT		 //mergeMode
                  , MPI_INT		 //sendMode
                  , MPI_INT		 //_algorithmSection
                  , MPI_CXX_BOOL		 //hasRefined
                  , MPI_CXX_BOOL		 //hasTriggeredRefinementForNextIteration
                  , MPI_CXX_BOOL		 //hasErased
                  , MPI_CXX_BOOL		 //hasTriggeredEraseForNextIteration
                  , MPI_CXX_BOOL		 //hasChangedVertexOrCellState
                  , MPI_CXX_BOOL		 //hasModifiedGridInPreviousIteration
                  , MPI_CXX_BOOL		 //isTraversalInverted
                  , MPI_CXX_BOOL		 //reduceStateAndCell
                  , MPI_CXX_BOOL		 //couldNotEraseDueToDecompositionFlag
                  , MPI_CXX_BOOL		 //subWorkerIsInvolvedInJoinOrFork
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //maxRefinementLevelAllowed
                  , 1		 //mergeMode
                  , 1		 //sendMode
                  , 1		 //_algorithmSection
                  , 1		 //hasRefined
                  , 1		 //hasTriggeredRefinementForNextIteration
                  , 1		 //hasErased
                  , 1		 //hasTriggeredEraseForNextIteration
                  , 1		 //hasChangedVertexOrCellState
                  , 1		 //hasModifiedGridInPreviousIteration
                  , 1		 //isTraversalInverted
                  , 1		 //reduceStateAndCell
                  , 1		 //couldNotEraseDueToDecompositionFlag
                  , 1		 //subWorkerIsInvolvedInJoinOrFork
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reduceStateAndCell))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reduceStateAndCell))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._couldNotEraseDueToDecompositionFlag))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._couldNotEraseDueToDecompositionFlag))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._subWorkerIsInvolvedInJoinOrFork))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._subWorkerIsInvolvedInJoinOrFork))), 		&disp[13] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(State)), i, disp[i], Attributes, sizeof(State));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]))), 		&disp[14] );
               disp[14] -= base;
               disp[14] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &State::Datatype );
               MPI_Type_commit( &State::Datatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::Datatype);
               MPI_Type_commit( &State::Datatype );
               #endif
               
            }
            {
               State dummyState[2];
               
               #ifdef MPI2
               const int Attributes = 14;
               #else
               const int Attributes = 15;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //maxRefinementLevelAllowed
                  , MPI_INT		 //mergeMode
                  , MPI_INT		 //sendMode
                  , MPI_INT		 //_algorithmSection
                  , MPI_CXX_BOOL		 //hasRefined
                  , MPI_CXX_BOOL		 //hasTriggeredRefinementForNextIteration
                  , MPI_CXX_BOOL		 //hasErased
                  , MPI_CXX_BOOL		 //hasTriggeredEraseForNextIteration
                  , MPI_CXX_BOOL		 //hasChangedVertexOrCellState
                  , MPI_CXX_BOOL		 //hasModifiedGridInPreviousIteration
                  , MPI_CXX_BOOL		 //isTraversalInverted
                  , MPI_CXX_BOOL		 //reduceStateAndCell
                  , MPI_CXX_BOOL		 //couldNotEraseDueToDecompositionFlag
                  , MPI_CXX_BOOL		 //subWorkerIsInvolvedInJoinOrFork
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //maxRefinementLevelAllowed
                  , 1		 //mergeMode
                  , 1		 //sendMode
                  , 1		 //_algorithmSection
                  , 1		 //hasRefined
                  , 1		 //hasTriggeredRefinementForNextIteration
                  , 1		 //hasErased
                  , 1		 //hasTriggeredEraseForNextIteration
                  , 1		 //hasChangedVertexOrCellState
                  , 1		 //hasModifiedGridInPreviousIteration
                  , 1		 //isTraversalInverted
                  , 1		 //reduceStateAndCell
                  , 1		 //couldNotEraseDueToDecompositionFlag
                  , 1		 //subWorkerIsInvolvedInJoinOrFork
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._sendMode))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasRefined))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredRefinementForNextIteration))), 		&disp[5] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[6] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasErased))), 		&disp[6] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[7] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasTriggeredEraseForNextIteration))), 		&disp[7] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[8] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasChangedVertexOrCellState))), 		&disp[8] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[9] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._hasModifiedGridInPreviousIteration))), 		&disp[9] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[10] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._isTraversalInverted))), 		&disp[10] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reduceStateAndCell))), 		&disp[11] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._reduceStateAndCell))), 		&disp[11] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._couldNotEraseDueToDecompositionFlag))), 		&disp[12] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._couldNotEraseDueToDecompositionFlag))), 		&disp[12] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._subWorkerIsInvolvedInJoinOrFork))), 		&disp[13] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[0]._persistentRecords._subWorkerIsInvolvedInJoinOrFork))), 		&disp[13] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(State)), i, disp[i], Attributes, sizeof(State));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyState[1]))), 		&disp[14] );
               disp[14] -= base;
               disp[14] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &State::FullDatatype );
               MPI_Type_commit( &State::FullDatatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &State::FullDatatype);
               MPI_Type_commit( &State::FullDatatype );
               #endif
               
            }
            
         }
         
         
         void exahype::records::State::shutdownDatatype() {
            MPI_Type_free( &State::Datatype );
            MPI_Type_free( &State::FullDatatype );
            
         }
         
         void exahype::records::State::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
            _senderDestinationRank = destination;
            
            if (communicateSleep<0) {
            
               const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
               if  (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "was not able to send message exahype::records::State "
                  << toString()
                  << " to node " << destination
                  << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "send(int)",msg.str() );
               }
               
            }
            else {
            
               MPI_Request* sendRequestHandle = new MPI_Request();
               MPI_Status   status;
               int          flag = 0;
               int          result;
               
               clock_t      timeOutWarning   = -1;
               clock_t      timeOutShutdown  = -1;
               bool         triggeredTimeoutWarning = false;
               
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  result = MPI_Isend(
                     this, 1, Datatype, destination,
                     tag, tarch::parallel::Node::getInstance().getCommunicator(),
                     sendRequestHandle
                  );
                  
               }
               else {
                  result = MPI_Isend(
                     this, 1, FullDatatype, destination,
                     tag, tarch::parallel::Node::getInstance().getCommunicator(),
                     sendRequestHandle
                  );
                  
               }
               if  (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "was not able to send message exahype::records::State "
                  << toString()
                  << " to node " << destination
                  << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "send(int)",msg.str() );
               }
               result = MPI_Test( sendRequestHandle, &flag, &status );
               while (!flag) {
                  if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
                  if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
                  result = MPI_Test( sendRequestHandle, &flag, &status );
                  if (result!=MPI_SUCCESS) {
                     std::ostringstream msg;
                     msg << "testing for finished send task for exahype::records::State "
                     << toString()
                     << " sent to node " << destination
                     << " failed: " << tarch::parallel::MPIReturnValueToString(result);
                     _log.error("send(int)", msg.str() );
                  }
                  
                  // deadlock aspect
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                     (clock()>timeOutWarning) &&
                     (!triggeredTimeoutWarning)
                  ) {
                     tarch::parallel::Node::getInstance().writeTimeOutWarning(
                     "exahype::records::State",
                     "send(int)", destination,tag,1
                     );
                     triggeredTimeoutWarning = true;
                  }
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                     (clock()>timeOutShutdown)
                  ) {
                     tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                     "exahype::records::State",
                     "send(int)", destination,tag,1
                     );
                  }
                  
               tarch::parallel::Node::getInstance().receiveDanglingMessages();
               usleep(communicateSleep);
               }
               
               delete sendRequestHandle;
               #ifdef Debug
               _log.debug("send(int,int)", "sent " + toString() );
               #endif
               
            }
            
         }
         
         
         
         void exahype::records::State::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
            if (communicateSleep<0) {
            
               MPI_Status  status;
               const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
               _senderDestinationRank = status.MPI_SOURCE;
               if ( result != MPI_SUCCESS ) {
                  std::ostringstream msg;
                  msg << "failed to start to receive exahype::records::State from node "
                  << source << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "receive(int)", msg.str() );
               }
               
            }
            else {
            
               MPI_Request* sendRequestHandle = new MPI_Request();
               MPI_Status   status;
               int          flag = 0;
               int          result;
               
               clock_t      timeOutWarning   = -1;
               clock_t      timeOutShutdown  = -1;
               bool         triggeredTimeoutWarning = false;
               
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  result = MPI_Irecv(
                     this, 1, Datatype, source, tag,
                     tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
                  );
                  
               }
               else {
                  result = MPI_Irecv(
                     this, 1, FullDatatype, source, tag,
                     tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
                  );
                  
               }
               if ( result != MPI_SUCCESS ) {
                  std::ostringstream msg;
                  msg << "failed to start to receive exahype::records::State from node "
                  << source << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "receive(int)", msg.str() );
               }
               
               result = MPI_Test( sendRequestHandle, &flag, &status );
               while (!flag) {
                  if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
                  if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
                  result = MPI_Test( sendRequestHandle, &flag, &status );
                  if (result!=MPI_SUCCESS) {
                     std::ostringstream msg;
                     msg << "testing for finished receive task for exahype::records::State failed: "
                     << tarch::parallel::MPIReturnValueToString(result);
                     _log.error("receive(int)", msg.str() );
                  }
                  
                  // deadlock aspect
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                     (clock()>timeOutWarning) &&
                     (!triggeredTimeoutWarning)
                  ) {
                     tarch::parallel::Node::getInstance().writeTimeOutWarning(
                     "exahype::records::State",
                     "receive(int)", source,tag,1
                     );
                     triggeredTimeoutWarning = true;
                  }
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                     (clock()>timeOutShutdown)
                  ) {
                     tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                     "exahype::records::State",
                     "receive(int)", source,tag,1
                     );
                  }
                  tarch::parallel::Node::getInstance().receiveDanglingMessages();
                  usleep(communicateSleep);
                  
               }
               
               delete sendRequestHandle;
               
               _senderDestinationRank = status.MPI_SOURCE;
               #ifdef Debug
               _log.debug("receive(int,int)", "received " + toString() ); 
               #endif
               
            }
            
         }
         
         
         
         bool exahype::records::State::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
            MPI_Status status;
            int  flag        = 0;
            MPI_Iprobe(
               MPI_ANY_SOURCE, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
            );
            if (flag) {
               int  messageCounter;
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  MPI_Get_count(&status, Datatype, &messageCounter);
               }
               else {
                  MPI_Get_count(&status, FullDatatype, &messageCounter);
               }
               return messageCounter > 0;
            }
            else return false;
            
         }
         
         int exahype::records::State::getSenderRank() const {
            assertion( _senderDestinationRank!=-1 );
            return _senderDestinationRank;
            
         }
      #endif
      
      
      exahype::records::StatePacked::PersistentRecords::PersistentRecords() {
         if ((9 >= (8 * sizeof(short int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((9 < (8 * sizeof(short int))));
         
      }
      
      
      exahype::records::StatePacked::PersistentRecords::PersistentRecords(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
      _maxRefinementLevelAllowed(maxRefinementLevelAllowed),
      _mergeMode(mergeMode),
      _sendMode(sendMode),
      _algorithmSection(algorithmSection),
      _isTraversalInverted(isTraversalInverted) {
         setHasRefined(hasRefined);
         setHasTriggeredRefinementForNextIteration(hasTriggeredRefinementForNextIteration);
         setHasErased(hasErased);
         setHasTriggeredEraseForNextIteration(hasTriggeredEraseForNextIteration);
         setHasChangedVertexOrCellState(hasChangedVertexOrCellState);
         setHasModifiedGridInPreviousIteration(hasModifiedGridInPreviousIteration);
         setReduceStateAndCell(reduceStateAndCell);
         setCouldNotEraseDueToDecompositionFlag(couldNotEraseDueToDecompositionFlag);
         setSubWorkerIsInvolvedInJoinOrFork(subWorkerIsInvolvedInJoinOrFork);
         if ((9 >= (8 * sizeof(short int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((9 < (8 * sizeof(short int))));
         
      }
      
      exahype::records::StatePacked::StatePacked() {
         if ((9 >= (8 * sizeof(short int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((9 < (8 * sizeof(short int))));
         
      }
      
      
      exahype::records::StatePacked::StatePacked(const PersistentRecords& persistentRecords):
      _persistentRecords(persistentRecords._maxRefinementLevelAllowed, persistentRecords._mergeMode, persistentRecords._sendMode, persistentRecords._algorithmSection, persistentRecords.getHasRefined(), persistentRecords.getHasTriggeredRefinementForNextIteration(), persistentRecords.getHasErased(), persistentRecords.getHasTriggeredEraseForNextIteration(), persistentRecords.getHasChangedVertexOrCellState(), persistentRecords.getHasModifiedGridInPreviousIteration(), persistentRecords._isTraversalInverted, persistentRecords.getReduceStateAndCell(), persistentRecords.getCouldNotEraseDueToDecompositionFlag(), persistentRecords.getSubWorkerIsInvolvedInJoinOrFork()) {
         if ((9 >= (8 * sizeof(short int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((9 < (8 * sizeof(short int))));
         
      }
      
      
      exahype::records::StatePacked::StatePacked(const int& maxRefinementLevelAllowed, const MergeMode& mergeMode, const SendMode& sendMode, const AlgorithmSection& algorithmSection, const bool& hasRefined, const bool& hasTriggeredRefinementForNextIteration, const bool& hasErased, const bool& hasTriggeredEraseForNextIteration, const bool& hasChangedVertexOrCellState, const bool& hasModifiedGridInPreviousIteration, const bool& isTraversalInverted, const bool& reduceStateAndCell, const bool& couldNotEraseDueToDecompositionFlag, const bool& subWorkerIsInvolvedInJoinOrFork):
      _persistentRecords(maxRefinementLevelAllowed, mergeMode, sendMode, algorithmSection, hasRefined, hasTriggeredRefinementForNextIteration, hasErased, hasTriggeredEraseForNextIteration, hasChangedVertexOrCellState, hasModifiedGridInPreviousIteration, isTraversalInverted, reduceStateAndCell, couldNotEraseDueToDecompositionFlag, subWorkerIsInvolvedInJoinOrFork) {
         if ((9 >= (8 * sizeof(short int)))) {
            std::cerr << "Packed-Type in " << __FILE__ << " too small. Either use bigger data type or append " << std::endl << std::endl;
            std::cerr << "  Packed-Type: short int hint-size no-of-bits;  " << std::endl << std::endl;
            std::cerr << "to your data type spec to guide DaStGen how many bits (no-of-bits) a data type has on your machine. DaStGen then can split up the bitfields into several attributes. " << std::endl; 
         }
         assertion((9 < (8 * sizeof(short int))));
         
      }
      
      
      exahype::records::StatePacked::~StatePacked() { }
      
      std::string exahype::records::StatePacked::toString(const MergeMode& param) {
         return exahype::records::State::toString(param);
      }
      
      std::string exahype::records::StatePacked::getMergeModeMapping() {
         return exahype::records::State::getMergeModeMapping();
      }
      
      std::string exahype::records::StatePacked::toString(const SendMode& param) {
         return exahype::records::State::toString(param);
      }
      
      std::string exahype::records::StatePacked::getSendModeMapping() {
         return exahype::records::State::getSendModeMapping();
      }
      
      std::string exahype::records::StatePacked::toString(const AlgorithmSection& param) {
         return exahype::records::State::toString(param);
      }
      
      std::string exahype::records::StatePacked::getAlgorithmSectionMapping() {
         return exahype::records::State::getAlgorithmSectionMapping();
      }
      
      
      
      std::string exahype::records::StatePacked::toString() const {
         std::ostringstream stringstr;
         toString(stringstr);
         return stringstr.str();
      }
      
      void exahype::records::StatePacked::toString (std::ostream& out) const {
         out << "("; 
         out << "maxRefinementLevelAllowed:" << getMaxRefinementLevelAllowed();
         out << ",";
         out << "mergeMode:" << toString(getMergeMode());
         out << ",";
         out << "sendMode:" << toString(getSendMode());
         out << ",";
         out << "_algorithmSection:" << toString(getAlgorithmSection());
         out << ",";
         out << "hasRefined:" << getHasRefined();
         out << ",";
         out << "hasTriggeredRefinementForNextIteration:" << getHasTriggeredRefinementForNextIteration();
         out << ",";
         out << "hasErased:" << getHasErased();
         out << ",";
         out << "hasTriggeredEraseForNextIteration:" << getHasTriggeredEraseForNextIteration();
         out << ",";
         out << "hasChangedVertexOrCellState:" << getHasChangedVertexOrCellState();
         out << ",";
         out << "hasModifiedGridInPreviousIteration:" << getHasModifiedGridInPreviousIteration();
         out << ",";
         out << "isTraversalInverted:" << getIsTraversalInverted();
         out << ",";
         out << "reduceStateAndCell:" << getReduceStateAndCell();
         out << ",";
         out << "couldNotEraseDueToDecompositionFlag:" << getCouldNotEraseDueToDecompositionFlag();
         out << ",";
         out << "subWorkerIsInvolvedInJoinOrFork:" << getSubWorkerIsInvolvedInJoinOrFork();
         out <<  ")";
      }
      
      
      exahype::records::StatePacked::PersistentRecords exahype::records::StatePacked::getPersistentRecords() const {
         return _persistentRecords;
      }
      
      exahype::records::State exahype::records::StatePacked::convert() const{
         return State(
            getMaxRefinementLevelAllowed(),
            getMergeMode(),
            getSendMode(),
            getAlgorithmSection(),
            getHasRefined(),
            getHasTriggeredRefinementForNextIteration(),
            getHasErased(),
            getHasTriggeredEraseForNextIteration(),
            getHasChangedVertexOrCellState(),
            getHasModifiedGridInPreviousIteration(),
            getIsTraversalInverted(),
            getReduceStateAndCell(),
            getCouldNotEraseDueToDecompositionFlag(),
            getSubWorkerIsInvolvedInJoinOrFork()
         );
      }
      
      #ifdef Parallel
         tarch::logging::Log exahype::records::StatePacked::_log( "exahype::records::StatePacked" );
         
         MPI_Datatype exahype::records::StatePacked::Datatype = 0;
         MPI_Datatype exahype::records::StatePacked::FullDatatype = 0;
         
         
         void exahype::records::StatePacked::initDatatype() {
            {
               StatePacked dummyStatePacked[2];
               
               #ifdef MPI2
               const int Attributes = 6;
               #else
               const int Attributes = 7;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //maxRefinementLevelAllowed
                  , MPI_INT		 //mergeMode
                  , MPI_INT		 //sendMode
                  , MPI_INT		 //_algorithmSection
                  , MPI_CXX_BOOL		 //isTraversalInverted
                  , MPI_SHORT		 //_packedRecords0
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //maxRefinementLevelAllowed
                  , 1		 //mergeMode
                  , 1		 //sendMode
                  , 1		 //_algorithmSection
                  , 1		 //isTraversalInverted
                  , 1		 //_packedRecords0
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[5] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(StatePacked)), i, disp[i], Attributes, sizeof(StatePacked));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]))), 		&disp[6] );
               disp[6] -= base;
               disp[6] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &StatePacked::Datatype );
               MPI_Type_commit( &StatePacked::Datatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::Datatype);
               MPI_Type_commit( &StatePacked::Datatype );
               #endif
               
            }
            {
               StatePacked dummyStatePacked[2];
               
               #ifdef MPI2
               const int Attributes = 6;
               #else
               const int Attributes = 7;
               #endif
               MPI_Datatype subtypes[Attributes] = {
                    MPI_INT		 //maxRefinementLevelAllowed
                  , MPI_INT		 //mergeMode
                  , MPI_INT		 //sendMode
                  , MPI_INT		 //_algorithmSection
                  , MPI_CXX_BOOL		 //isTraversalInverted
                  , MPI_SHORT		 //_packedRecords0
                  #ifndef MPI2
                  , MPI_UB
                  #endif
                  
               };
               
               int blocklen[Attributes] = {
                    1		 //maxRefinementLevelAllowed
                  , 1		 //mergeMode
                  , 1		 //sendMode
                  , 1		 //_algorithmSection
                  , 1		 //isTraversalInverted
                  , 1		 //_packedRecords0
                  #ifndef MPI2
                  , 1
                  #endif
                  
               };
               
               MPI_Aint  disp[Attributes];
               MPI_Aint  base;
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked))), &base);
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._maxRefinementLevelAllowed))), 		&disp[0] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._mergeMode))), 		&disp[1] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._sendMode))), 		&disp[2] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._algorithmSection))), 		&disp[3] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[4] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._isTraversalInverted))), 		&disp[4] );
               #endif
               #ifdef MPI2
               MPI_Get_address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[5] );
               #else
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[0]._persistentRecords._packedRecords0))), 		&disp[5] );
               #endif
               #ifdef MPI2
               for (int i=1; i<Attributes; i++) {
               #else
               for (int i=1; i<Attributes-1; i++) {
               #endif
                  assertion1( disp[i] > disp[i-1], i );
               }
               #ifdef MPI2
               for (int i=0; i<Attributes; i++) {
               #else
               for (int i=0; i<Attributes-1; i++) {
               #endif
                  disp[i] = disp[i] - base; // should be MPI_Aint_diff(disp[i], base); but this is not supported by most MPI-2 implementations
                  assertion4(disp[i]<static_cast<int>(sizeof(StatePacked)), i, disp[i], Attributes, sizeof(StatePacked));
               }
               #ifndef MPI2
               MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyStatePacked[1]))), 		&disp[6] );
               disp[6] -= base;
               disp[6] += disp[0];
               #endif
               #ifdef MPI2
               MPI_Datatype tmpType; 
               MPI_Aint lowerBound, typeExtent; 
               MPI_Type_create_struct( Attributes, blocklen, disp, subtypes, &tmpType );
               MPI_Type_get_extent( tmpType, &lowerBound, &typeExtent );
               MPI_Type_create_resized( tmpType, lowerBound, typeExtent, &StatePacked::FullDatatype );
               MPI_Type_commit( &StatePacked::FullDatatype );
               #else
               MPI_Type_struct( Attributes, blocklen, disp, subtypes, &StatePacked::FullDatatype);
               MPI_Type_commit( &StatePacked::FullDatatype );
               #endif
               
            }
            
         }
         
         
         void exahype::records::StatePacked::shutdownDatatype() {
            MPI_Type_free( &StatePacked::Datatype );
            MPI_Type_free( &StatePacked::FullDatatype );
            
         }
         
         void exahype::records::StatePacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
            _senderDestinationRank = destination;
            
            if (communicateSleep<0) {
            
               const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
               if  (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "was not able to send message exahype::records::StatePacked "
                  << toString()
                  << " to node " << destination
                  << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "send(int)",msg.str() );
               }
               
            }
            else {
            
               MPI_Request* sendRequestHandle = new MPI_Request();
               MPI_Status   status;
               int          flag = 0;
               int          result;
               
               clock_t      timeOutWarning   = -1;
               clock_t      timeOutShutdown  = -1;
               bool         triggeredTimeoutWarning = false;
               
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  result = MPI_Isend(
                     this, 1, Datatype, destination,
                     tag, tarch::parallel::Node::getInstance().getCommunicator(),
                     sendRequestHandle
                  );
                  
               }
               else {
                  result = MPI_Isend(
                     this, 1, FullDatatype, destination,
                     tag, tarch::parallel::Node::getInstance().getCommunicator(),
                     sendRequestHandle
                  );
                  
               }
               if  (result!=MPI_SUCCESS) {
                  std::ostringstream msg;
                  msg << "was not able to send message exahype::records::StatePacked "
                  << toString()
                  << " to node " << destination
                  << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "send(int)",msg.str() );
               }
               result = MPI_Test( sendRequestHandle, &flag, &status );
               while (!flag) {
                  if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
                  if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
                  result = MPI_Test( sendRequestHandle, &flag, &status );
                  if (result!=MPI_SUCCESS) {
                     std::ostringstream msg;
                     msg << "testing for finished send task for exahype::records::StatePacked "
                     << toString()
                     << " sent to node " << destination
                     << " failed: " << tarch::parallel::MPIReturnValueToString(result);
                     _log.error("send(int)", msg.str() );
                  }
                  
                  // deadlock aspect
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                     (clock()>timeOutWarning) &&
                     (!triggeredTimeoutWarning)
                  ) {
                     tarch::parallel::Node::getInstance().writeTimeOutWarning(
                     "exahype::records::StatePacked",
                     "send(int)", destination,tag,1
                     );
                     triggeredTimeoutWarning = true;
                  }
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                     (clock()>timeOutShutdown)
                  ) {
                     tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                     "exahype::records::StatePacked",
                     "send(int)", destination,tag,1
                     );
                  }
                  
               tarch::parallel::Node::getInstance().receiveDanglingMessages();
               usleep(communicateSleep);
               }
               
               delete sendRequestHandle;
               #ifdef Debug
               _log.debug("send(int,int)", "sent " + toString() );
               #endif
               
            }
            
         }
         
         
         
         void exahype::records::StatePacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
            if (communicateSleep<0) {
            
               MPI_Status  status;
               const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
               _senderDestinationRank = status.MPI_SOURCE;
               if ( result != MPI_SUCCESS ) {
                  std::ostringstream msg;
                  msg << "failed to start to receive exahype::records::StatePacked from node "
                  << source << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "receive(int)", msg.str() );
               }
               
            }
            else {
            
               MPI_Request* sendRequestHandle = new MPI_Request();
               MPI_Status   status;
               int          flag = 0;
               int          result;
               
               clock_t      timeOutWarning   = -1;
               clock_t      timeOutShutdown  = -1;
               bool         triggeredTimeoutWarning = false;
               
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  result = MPI_Irecv(
                     this, 1, Datatype, source, tag,
                     tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
                  );
                  
               }
               else {
                  result = MPI_Irecv(
                     this, 1, FullDatatype, source, tag,
                     tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
                  );
                  
               }
               if ( result != MPI_SUCCESS ) {
                  std::ostringstream msg;
                  msg << "failed to start to receive exahype::records::StatePacked from node "
                  << source << ": " << tarch::parallel::MPIReturnValueToString(result);
                  _log.error( "receive(int)", msg.str() );
               }
               
               result = MPI_Test( sendRequestHandle, &flag, &status );
               while (!flag) {
                  if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
                  if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
                  result = MPI_Test( sendRequestHandle, &flag, &status );
                  if (result!=MPI_SUCCESS) {
                     std::ostringstream msg;
                     msg << "testing for finished receive task for exahype::records::StatePacked failed: "
                     << tarch::parallel::MPIReturnValueToString(result);
                     _log.error("receive(int)", msg.str() );
                  }
                  
                  // deadlock aspect
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
                     (clock()>timeOutWarning) &&
                     (!triggeredTimeoutWarning)
                  ) {
                     tarch::parallel::Node::getInstance().writeTimeOutWarning(
                     "exahype::records::StatePacked",
                     "receive(int)", source,tag,1
                     );
                     triggeredTimeoutWarning = true;
                  }
                  if (
                     tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
                     (clock()>timeOutShutdown)
                  ) {
                     tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
                     "exahype::records::StatePacked",
                     "receive(int)", source,tag,1
                     );
                  }
                  tarch::parallel::Node::getInstance().receiveDanglingMessages();
                  usleep(communicateSleep);
                  
               }
               
               delete sendRequestHandle;
               
               _senderDestinationRank = status.MPI_SOURCE;
               #ifdef Debug
               _log.debug("receive(int,int)", "received " + toString() ); 
               #endif
               
            }
            
         }
         
         
         
         bool exahype::records::StatePacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
            MPI_Status status;
            int  flag        = 0;
            MPI_Iprobe(
               MPI_ANY_SOURCE, tag,
               tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
            );
            if (flag) {
               int  messageCounter;
               if (exchangeOnlyAttributesMarkedWithParallelise) {
                  MPI_Get_count(&status, Datatype, &messageCounter);
               }
               else {
                  MPI_Get_count(&status, FullDatatype, &messageCounter);
               }
               return messageCounter > 0;
            }
            else return false;
            
         }
         
         int exahype::records::StatePacked::getSenderRank() const {
            assertion( _senderDestinationRank!=-1 );
            return _senderDestinationRank;
            
         }
      #endif
      
      
      
   
#endif


