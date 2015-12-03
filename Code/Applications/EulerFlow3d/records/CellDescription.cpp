#include "EulerFlow3d/records/CellDescription.h"

exahype::records::CellDescription::PersistentRecords::PersistentRecords() {
   
}


exahype::records::CellDescription::PersistentRecords::PersistentRecords(const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& spaceTimePredictor, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& spaceTimeVolumeFlux, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& solution, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& update, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& predictor, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& volumeFlux, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& extrapolatedPredictor, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& fluctuation, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size):
_spaceTimePredictor(spaceTimePredictor),
_spaceTimeVolumeFlux(spaceTimeVolumeFlux),
_solution(solution),
_update(update),
_predictor(predictor),
_volumeFlux(volumeFlux),
_extrapolatedPredictor(extrapolatedPredictor),
_fluctuation(fluctuation),
_level(level),
_offset(offset),
_size(size) {
   
}

exahype::records::CellDescription::CellDescription() {
   
}


exahype::records::CellDescription::CellDescription(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._spaceTimePredictor, persistentRecords._spaceTimeVolumeFlux, persistentRecords._solution, persistentRecords._update, persistentRecords._predictor, persistentRecords._volumeFlux, persistentRecords._extrapolatedPredictor, persistentRecords._fluctuation, persistentRecords._level, persistentRecords._offset, persistentRecords._size) {
   
}


exahype::records::CellDescription::CellDescription(const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& spaceTimePredictor, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& spaceTimeVolumeFlux, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& solution, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& update, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& predictor, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& volumeFlux, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& extrapolatedPredictor, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& fluctuation, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size):
_persistentRecords(spaceTimePredictor, spaceTimeVolumeFlux, solution, update, predictor, volumeFlux, extrapolatedPredictor, fluctuation, level, offset, size) {
   
}


exahype::records::CellDescription::~CellDescription() { }



std::string exahype::records::CellDescription::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void exahype::records::CellDescription::toString (std::ostream& out) const {
   out << "("; 
   out << "spaceTimePredictor:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getSpaceTimePredictor(i) << ",";
   }
   out << getSpaceTimePredictor(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "spaceTimeVolumeFlux:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getSpaceTimeVolumeFlux(i) << ",";
   }
   out << getSpaceTimeVolumeFlux(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "solution:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getSolution(i) << ",";
   }
   out << getSolution(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "update:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getUpdate(i) << ",";
   }
   out << getUpdate(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "predictor:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getPredictor(i) << ",";
   }
   out << getPredictor(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "volumeFlux:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getVolumeFlux(i) << ",";
   }
   out << getVolumeFlux(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "extrapolatedPredictor:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getExtrapolatedPredictor(i) << ",";
   }
   out << getExtrapolatedPredictor(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "fluctuation:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getFluctuation(i) << ",";
   }
   out << getFluctuation(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "level:" << getLevel();
   out << ",";
   out << "offset:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getOffset(i) << ",";
   }
   out << getOffset(DIMENSIONS-1) << "]";
   out << ",";
   out << "size:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getSize(i) << ",";
   }
   out << getSize(DIMENSIONS-1) << "]";
   out <<  ")";
}


exahype::records::CellDescription::PersistentRecords exahype::records::CellDescription::getPersistentRecords() const {
   return _persistentRecords;
}

exahype::records::CellDescriptionPacked exahype::records::CellDescription::convert() const{
   return CellDescriptionPacked(
      getSpaceTimePredictor(),
      getSpaceTimeVolumeFlux(),
      getSolution(),
      getUpdate(),
      getPredictor(),
      getVolumeFlux(),
      getExtrapolatedPredictor(),
      getFluctuation(),
      getLevel(),
      getOffset(),
      getSize()
   );
}

#ifdef Parallel
   tarch::logging::Log exahype::records::CellDescription::_log( "exahype::records::CellDescription" );
   
   MPI_Datatype exahype::records::CellDescription::Datatype = 0;
   MPI_Datatype exahype::records::CellDescription::FullDatatype = 0;
   
   
   void exahype::records::CellDescription::initDatatype() {
      {
         CellDescription dummyCellDescription[2];
         
         const int Attributes = 12;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //spaceTimePredictor
            MPI_INT,		 //spaceTimeVolumeFlux
            MPI_INT,		 //solution
            MPI_INT,		 //update
            MPI_INT,		 //predictor
            MPI_INT,		 //volumeFlux
            MPI_INT,		 //extrapolatedPredictor
            MPI_INT,		 //fluctuation
            MPI_INT,		 //level
            MPI_DOUBLE,		 //offset
            MPI_DOUBLE,		 //size
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            EXAHYPE_PATCH_SIZE_TOTAL,		 //spaceTimePredictor
            EXAHYPE_PATCH_SIZE_TOTAL,		 //spaceTimeVolumeFlux
            EXAHYPE_PATCH_SIZE_TOTAL,		 //solution
            EXAHYPE_PATCH_SIZE_TOTAL,		 //update
            EXAHYPE_PATCH_SIZE_TOTAL,		 //predictor
            EXAHYPE_PATCH_SIZE_TOTAL,		 //volumeFlux
            EXAHYPE_PATCH_SIZE_TOTAL,		 //extrapolatedPredictor
            EXAHYPE_PATCH_SIZE_TOTAL,		 //fluctuation
            1,		 //level
            DIMENSIONS,		 //offset
            DIMENSIONS,		 //size
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._spaceTimePredictor[0]))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._spaceTimeVolumeFlux[0]))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._solution[0]))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._update[0]))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._predictor[0]))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._volumeFlux[0]))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._extrapolatedPredictor[0]))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._fluctuation[0]))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._level))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._offset[0]))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._size[0]))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&dummyCellDescription[1]._persistentRecords._spaceTimePredictor[0])), 		&disp[11] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellDescription::Datatype );
         MPI_Type_commit( &CellDescription::Datatype );
         
      }
      {
         CellDescription dummyCellDescription[2];
         
         const int Attributes = 12;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //spaceTimePredictor
            MPI_INT,		 //spaceTimeVolumeFlux
            MPI_INT,		 //solution
            MPI_INT,		 //update
            MPI_INT,		 //predictor
            MPI_INT,		 //volumeFlux
            MPI_INT,		 //extrapolatedPredictor
            MPI_INT,		 //fluctuation
            MPI_INT,		 //level
            MPI_DOUBLE,		 //offset
            MPI_DOUBLE,		 //size
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            EXAHYPE_PATCH_SIZE_TOTAL,		 //spaceTimePredictor
            EXAHYPE_PATCH_SIZE_TOTAL,		 //spaceTimeVolumeFlux
            EXAHYPE_PATCH_SIZE_TOTAL,		 //solution
            EXAHYPE_PATCH_SIZE_TOTAL,		 //update
            EXAHYPE_PATCH_SIZE_TOTAL,		 //predictor
            EXAHYPE_PATCH_SIZE_TOTAL,		 //volumeFlux
            EXAHYPE_PATCH_SIZE_TOTAL,		 //extrapolatedPredictor
            EXAHYPE_PATCH_SIZE_TOTAL,		 //fluctuation
            1,		 //level
            DIMENSIONS,		 //offset
            DIMENSIONS,		 //size
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._spaceTimePredictor[0]))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._spaceTimeVolumeFlux[0]))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._solution[0]))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._update[0]))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._predictor[0]))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._volumeFlux[0]))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._extrapolatedPredictor[0]))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._fluctuation[0]))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._level))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._offset[0]))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescription[0]._persistentRecords._size[0]))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&dummyCellDescription[1]._persistentRecords._spaceTimePredictor[0])), 		&disp[11] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellDescription::FullDatatype );
         MPI_Type_commit( &CellDescription::FullDatatype );
         
      }
      
   }
   
   
   void exahype::records::CellDescription::shutdownDatatype() {
      MPI_Type_free( &CellDescription::Datatype );
      MPI_Type_free( &CellDescription::FullDatatype );
      
   }
   
   void exahype::records::CellDescription::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      if (communicateSleep<0) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message exahype::records::CellDescription "
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
            msg << "was not able to send message exahype::records::CellDescription "
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
               msg << "testing for finished send task for exahype::records::CellDescription "
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
               "exahype::records::CellDescription",
               "send(int)", destination,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::CellDescription",
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
   
   
   
   void exahype::records::CellDescription::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      if (communicateSleep<0) {
      
         MPI_Status  status;
         const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive exahype::records::CellDescription from node "
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
            msg << "failed to start to receive exahype::records::CellDescription from node "
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
               msg << "testing for finished receive task for exahype::records::CellDescription failed: "
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
               "exahype::records::CellDescription",
               "receive(int)", source,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::CellDescription",
               "receive(int)", source,tag,1
               );
            }
            tarch::parallel::Node::getInstance().receiveDanglingMessages();
            usleep(communicateSleep);
            
         }
         
         delete sendRequestHandle;
         
         #ifdef Debug
         _log.debug("receive(int,int)", "received " + toString() ); 
         #endif
         
      }
      
   }
   
   
   
   bool exahype::records::CellDescription::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
   
   
#endif


exahype::records::CellDescriptionPacked::PersistentRecords::PersistentRecords() {
   
}


exahype::records::CellDescriptionPacked::PersistentRecords::PersistentRecords(const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& spaceTimePredictor, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& spaceTimeVolumeFlux, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& solution, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& update, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& predictor, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& volumeFlux, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& extrapolatedPredictor, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& fluctuation, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size):
_spaceTimePredictor(spaceTimePredictor),
_spaceTimeVolumeFlux(spaceTimeVolumeFlux),
_solution(solution),
_update(update),
_predictor(predictor),
_volumeFlux(volumeFlux),
_extrapolatedPredictor(extrapolatedPredictor),
_fluctuation(fluctuation),
_level(level),
_offset(offset),
_size(size) {
   
}

exahype::records::CellDescriptionPacked::CellDescriptionPacked() {
   
}


exahype::records::CellDescriptionPacked::CellDescriptionPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._spaceTimePredictor, persistentRecords._spaceTimeVolumeFlux, persistentRecords._solution, persistentRecords._update, persistentRecords._predictor, persistentRecords._volumeFlux, persistentRecords._extrapolatedPredictor, persistentRecords._fluctuation, persistentRecords._level, persistentRecords._offset, persistentRecords._size) {
   
}


exahype::records::CellDescriptionPacked::CellDescriptionPacked(const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& spaceTimePredictor, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& spaceTimeVolumeFlux, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& solution, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& update, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& predictor, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& volumeFlux, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& extrapolatedPredictor, const tarch::la::Vector<EXAHYPE_PATCH_SIZE_TOTAL,int>& fluctuation, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size):
_persistentRecords(spaceTimePredictor, spaceTimeVolumeFlux, solution, update, predictor, volumeFlux, extrapolatedPredictor, fluctuation, level, offset, size) {
   
}


exahype::records::CellDescriptionPacked::~CellDescriptionPacked() { }



std::string exahype::records::CellDescriptionPacked::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void exahype::records::CellDescriptionPacked::toString (std::ostream& out) const {
   out << "("; 
   out << "spaceTimePredictor:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getSpaceTimePredictor(i) << ",";
   }
   out << getSpaceTimePredictor(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "spaceTimeVolumeFlux:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getSpaceTimeVolumeFlux(i) << ",";
   }
   out << getSpaceTimeVolumeFlux(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "solution:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getSolution(i) << ",";
   }
   out << getSolution(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "update:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getUpdate(i) << ",";
   }
   out << getUpdate(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "predictor:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getPredictor(i) << ",";
   }
   out << getPredictor(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "volumeFlux:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getVolumeFlux(i) << ",";
   }
   out << getVolumeFlux(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "extrapolatedPredictor:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getExtrapolatedPredictor(i) << ",";
   }
   out << getExtrapolatedPredictor(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "fluctuation:[";
   for (int i = 0; i < EXAHYPE_PATCH_SIZE_TOTAL-1; i++) {
      out << getFluctuation(i) << ",";
   }
   out << getFluctuation(EXAHYPE_PATCH_SIZE_TOTAL-1) << "]";
   out << ",";
   out << "level:" << getLevel();
   out << ",";
   out << "offset:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getOffset(i) << ",";
   }
   out << getOffset(DIMENSIONS-1) << "]";
   out << ",";
   out << "size:[";
   for (int i = 0; i < DIMENSIONS-1; i++) {
      out << getSize(i) << ",";
   }
   out << getSize(DIMENSIONS-1) << "]";
   out <<  ")";
}


exahype::records::CellDescriptionPacked::PersistentRecords exahype::records::CellDescriptionPacked::getPersistentRecords() const {
   return _persistentRecords;
}

exahype::records::CellDescription exahype::records::CellDescriptionPacked::convert() const{
   return CellDescription(
      getSpaceTimePredictor(),
      getSpaceTimeVolumeFlux(),
      getSolution(),
      getUpdate(),
      getPredictor(),
      getVolumeFlux(),
      getExtrapolatedPredictor(),
      getFluctuation(),
      getLevel(),
      getOffset(),
      getSize()
   );
}

#ifdef Parallel
   tarch::logging::Log exahype::records::CellDescriptionPacked::_log( "exahype::records::CellDescriptionPacked" );
   
   MPI_Datatype exahype::records::CellDescriptionPacked::Datatype = 0;
   MPI_Datatype exahype::records::CellDescriptionPacked::FullDatatype = 0;
   
   
   void exahype::records::CellDescriptionPacked::initDatatype() {
      {
         CellDescriptionPacked dummyCellDescriptionPacked[2];
         
         const int Attributes = 12;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //spaceTimePredictor
            MPI_INT,		 //spaceTimeVolumeFlux
            MPI_INT,		 //solution
            MPI_INT,		 //update
            MPI_INT,		 //predictor
            MPI_INT,		 //volumeFlux
            MPI_INT,		 //extrapolatedPredictor
            MPI_INT,		 //fluctuation
            MPI_INT,		 //level
            MPI_DOUBLE,		 //offset
            MPI_DOUBLE,		 //size
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            EXAHYPE_PATCH_SIZE_TOTAL,		 //spaceTimePredictor
            EXAHYPE_PATCH_SIZE_TOTAL,		 //spaceTimeVolumeFlux
            EXAHYPE_PATCH_SIZE_TOTAL,		 //solution
            EXAHYPE_PATCH_SIZE_TOTAL,		 //update
            EXAHYPE_PATCH_SIZE_TOTAL,		 //predictor
            EXAHYPE_PATCH_SIZE_TOTAL,		 //volumeFlux
            EXAHYPE_PATCH_SIZE_TOTAL,		 //extrapolatedPredictor
            EXAHYPE_PATCH_SIZE_TOTAL,		 //fluctuation
            1,		 //level
            DIMENSIONS,		 //offset
            DIMENSIONS,		 //size
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._spaceTimePredictor[0]))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._spaceTimeVolumeFlux[0]))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._solution[0]))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._update[0]))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._predictor[0]))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._volumeFlux[0]))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictor[0]))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._fluctuation[0]))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&dummyCellDescriptionPacked[1]._persistentRecords._spaceTimePredictor[0])), 		&disp[11] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellDescriptionPacked::Datatype );
         MPI_Type_commit( &CellDescriptionPacked::Datatype );
         
      }
      {
         CellDescriptionPacked dummyCellDescriptionPacked[2];
         
         const int Attributes = 12;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //spaceTimePredictor
            MPI_INT,		 //spaceTimeVolumeFlux
            MPI_INT,		 //solution
            MPI_INT,		 //update
            MPI_INT,		 //predictor
            MPI_INT,		 //volumeFlux
            MPI_INT,		 //extrapolatedPredictor
            MPI_INT,		 //fluctuation
            MPI_INT,		 //level
            MPI_DOUBLE,		 //offset
            MPI_DOUBLE,		 //size
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            EXAHYPE_PATCH_SIZE_TOTAL,		 //spaceTimePredictor
            EXAHYPE_PATCH_SIZE_TOTAL,		 //spaceTimeVolumeFlux
            EXAHYPE_PATCH_SIZE_TOTAL,		 //solution
            EXAHYPE_PATCH_SIZE_TOTAL,		 //update
            EXAHYPE_PATCH_SIZE_TOTAL,		 //predictor
            EXAHYPE_PATCH_SIZE_TOTAL,		 //volumeFlux
            EXAHYPE_PATCH_SIZE_TOTAL,		 //extrapolatedPredictor
            EXAHYPE_PATCH_SIZE_TOTAL,		 //fluctuation
            1,		 //level
            DIMENSIONS,		 //offset
            DIMENSIONS,		 //size
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._spaceTimePredictor[0]))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._spaceTimeVolumeFlux[0]))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._solution[0]))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._update[0]))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._predictor[0]))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._volumeFlux[0]))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._extrapolatedPredictor[0]))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._fluctuation[0]))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&dummyCellDescriptionPacked[1]._persistentRecords._spaceTimePredictor[0])), 		&disp[11] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &CellDescriptionPacked::FullDatatype );
         MPI_Type_commit( &CellDescriptionPacked::FullDatatype );
         
      }
      
   }
   
   
   void exahype::records::CellDescriptionPacked::shutdownDatatype() {
      MPI_Type_free( &CellDescriptionPacked::Datatype );
      MPI_Type_free( &CellDescriptionPacked::FullDatatype );
      
   }
   
   void exahype::records::CellDescriptionPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      if (communicateSleep<0) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message exahype::records::CellDescriptionPacked "
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
            msg << "was not able to send message exahype::records::CellDescriptionPacked "
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
               msg << "testing for finished send task for exahype::records::CellDescriptionPacked "
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
               "exahype::records::CellDescriptionPacked",
               "send(int)", destination,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::CellDescriptionPacked",
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
   
   
   
   void exahype::records::CellDescriptionPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      if (communicateSleep<0) {
      
         MPI_Status  status;
         const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive exahype::records::CellDescriptionPacked from node "
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
            msg << "failed to start to receive exahype::records::CellDescriptionPacked from node "
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
               msg << "testing for finished receive task for exahype::records::CellDescriptionPacked failed: "
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
               "exahype::records::CellDescriptionPacked",
               "receive(int)", source,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::CellDescriptionPacked",
               "receive(int)", source,tag,1
               );
            }
            tarch::parallel::Node::getInstance().receiveDanglingMessages();
            usleep(communicateSleep);
            
         }
         
         delete sendRequestHandle;
         
         #ifdef Debug
         _log.debug("receive(int,int)", "received " + toString() ); 
         #endif
         
      }
      
   }
   
   
   
   bool exahype::records::CellDescriptionPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
   
   
#endif



