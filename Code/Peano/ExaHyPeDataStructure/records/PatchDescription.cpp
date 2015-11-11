#include "ExaHyPeDataStructure/records/PatchDescription.h"

exahype::records::PatchDescription::PersistentRecords::PersistentRecords() {
   
}


exahype::records::PatchDescription::PersistentRecords::PersistentRecords(const int& spaceTimePredictorDof, const int& spaceTimeVolumeFluxesDof, const int& dof, const int& updateDof, const int& predictorDof, const int& volumeFluxesDof, const int& extrapolatedPredictorDof, const int& normalFluxesDof, const int& fluctuationsDof, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size):
_spaceTimePredictorDof(spaceTimePredictorDof),
_spaceTimeVolumeFluxesDof(spaceTimeVolumeFluxesDof),
_dof(dof),
_updateDof(updateDof),
_predictorDof(predictorDof),
_volumeFluxesDof(volumeFluxesDof),
_extrapolatedPredictorDof(extrapolatedPredictorDof),
_normalFluxesDof(normalFluxesDof),
_fluctuationsDof(fluctuationsDof),
_level(level),
_offset(offset),
_size(size) {
   
}

exahype::records::PatchDescription::PatchDescription() {
   
}


exahype::records::PatchDescription::PatchDescription(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._spaceTimePredictorDof, persistentRecords._spaceTimeVolumeFluxesDof, persistentRecords._dof, persistentRecords._updateDof, persistentRecords._predictorDof, persistentRecords._volumeFluxesDof, persistentRecords._extrapolatedPredictorDof, persistentRecords._normalFluxesDof, persistentRecords._fluctuationsDof, persistentRecords._level, persistentRecords._offset, persistentRecords._size) {
   
}


exahype::records::PatchDescription::PatchDescription(const int& spaceTimePredictorDof, const int& spaceTimeVolumeFluxesDof, const int& dof, const int& updateDof, const int& predictorDof, const int& volumeFluxesDof, const int& extrapolatedPredictorDof, const int& normalFluxesDof, const int& fluctuationsDof, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size):
_persistentRecords(spaceTimePredictorDof, spaceTimeVolumeFluxesDof, dof, updateDof, predictorDof, volumeFluxesDof, extrapolatedPredictorDof, normalFluxesDof, fluctuationsDof, level, offset, size) {
   
}


exahype::records::PatchDescription::~PatchDescription() { }



std::string exahype::records::PatchDescription::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void exahype::records::PatchDescription::toString (std::ostream& out) const {
   out << "("; 
   out << "spaceTimePredictorDof:" << getSpaceTimePredictorDof();
   out << ",";
   out << "spaceTimeVolumeFluxesDof:" << getSpaceTimeVolumeFluxesDof();
   out << ",";
   out << "dof:" << getDof();
   out << ",";
   out << "updateDof:" << getUpdateDof();
   out << ",";
   out << "predictorDof:" << getPredictorDof();
   out << ",";
   out << "volumeFluxesDof:" << getVolumeFluxesDof();
   out << ",";
   out << "extrapolatedPredictorDof:" << getExtrapolatedPredictorDof();
   out << ",";
   out << "normalFluxesDof:" << getNormalFluxesDof();
   out << ",";
   out << "fluctuationsDof:" << getFluctuationsDof();
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


exahype::records::PatchDescription::PersistentRecords exahype::records::PatchDescription::getPersistentRecords() const {
   return _persistentRecords;
}

exahype::records::PatchDescriptionPacked exahype::records::PatchDescription::convert() const{
   return PatchDescriptionPacked(
      getSpaceTimePredictorDof(),
      getSpaceTimeVolumeFluxesDof(),
      getDof(),
      getUpdateDof(),
      getPredictorDof(),
      getVolumeFluxesDof(),
      getExtrapolatedPredictorDof(),
      getNormalFluxesDof(),
      getFluctuationsDof(),
      getLevel(),
      getOffset(),
      getSize()
   );
}

#ifdef Parallel
   tarch::logging::Log exahype::records::PatchDescription::_log( "exahype::records::PatchDescription" );
   
   MPI_Datatype exahype::records::PatchDescription::Datatype = 0;
   MPI_Datatype exahype::records::PatchDescription::FullDatatype = 0;
   
   
   void exahype::records::PatchDescription::initDatatype() {
      {
         PatchDescription dummyPatchDescription[2];
         
         const int Attributes = 13;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //spaceTimePredictorDof
            MPI_INT,		 //spaceTimeVolumeFluxesDof
            MPI_INT,		 //dof
            MPI_INT,		 //updateDof
            MPI_INT,		 //predictorDof
            MPI_INT,		 //volumeFluxesDof
            MPI_INT,		 //extrapolatedPredictorDof
            MPI_INT,		 //normalFluxesDof
            MPI_INT,		 //fluctuationsDof
            MPI_INT,		 //level
            MPI_DOUBLE,		 //offset
            MPI_DOUBLE,		 //size
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //spaceTimePredictorDof
            1,		 //spaceTimeVolumeFluxesDof
            1,		 //dof
            1,		 //updateDof
            1,		 //predictorDof
            1,		 //volumeFluxesDof
            1,		 //extrapolatedPredictorDof
            1,		 //normalFluxesDof
            1,		 //fluctuationsDof
            1,		 //level
            DIMENSIONS,		 //offset
            DIMENSIONS,		 //size
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._spaceTimePredictorDof))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._spaceTimeVolumeFluxesDof))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._dof))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._updateDof))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._predictorDof))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._volumeFluxesDof))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._extrapolatedPredictorDof))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._normalFluxesDof))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._fluctuationsDof))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._level))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._offset[0]))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._size[0]))), 		&disp[11] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[1]._persistentRecords._spaceTimePredictorDof))), 		&disp[12] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &PatchDescription::Datatype );
         MPI_Type_commit( &PatchDescription::Datatype );
         
      }
      {
         PatchDescription dummyPatchDescription[2];
         
         const int Attributes = 13;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //spaceTimePredictorDof
            MPI_INT,		 //spaceTimeVolumeFluxesDof
            MPI_INT,		 //dof
            MPI_INT,		 //updateDof
            MPI_INT,		 //predictorDof
            MPI_INT,		 //volumeFluxesDof
            MPI_INT,		 //extrapolatedPredictorDof
            MPI_INT,		 //normalFluxesDof
            MPI_INT,		 //fluctuationsDof
            MPI_INT,		 //level
            MPI_DOUBLE,		 //offset
            MPI_DOUBLE,		 //size
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //spaceTimePredictorDof
            1,		 //spaceTimeVolumeFluxesDof
            1,		 //dof
            1,		 //updateDof
            1,		 //predictorDof
            1,		 //volumeFluxesDof
            1,		 //extrapolatedPredictorDof
            1,		 //normalFluxesDof
            1,		 //fluctuationsDof
            1,		 //level
            DIMENSIONS,		 //offset
            DIMENSIONS,		 //size
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._spaceTimePredictorDof))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._spaceTimeVolumeFluxesDof))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._dof))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._updateDof))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._predictorDof))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._volumeFluxesDof))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._extrapolatedPredictorDof))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._normalFluxesDof))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._fluctuationsDof))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._level))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._offset[0]))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[0]._persistentRecords._size[0]))), 		&disp[11] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescription[1]._persistentRecords._spaceTimePredictorDof))), 		&disp[12] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &PatchDescription::FullDatatype );
         MPI_Type_commit( &PatchDescription::FullDatatype );
         
      }
      
   }
   
   
   void exahype::records::PatchDescription::shutdownDatatype() {
      MPI_Type_free( &PatchDescription::Datatype );
      MPI_Type_free( &PatchDescription::FullDatatype );
      
   }
   
   void exahype::records::PatchDescription::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      if (communicateSleep<0) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message exahype::records::PatchDescription "
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
            msg << "was not able to send message exahype::records::PatchDescription "
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
               msg << "testing for finished send task for exahype::records::PatchDescription "
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
               "exahype::records::PatchDescription",
               "send(int)", destination,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::PatchDescription",
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
   
   
   
   void exahype::records::PatchDescription::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      if (communicateSleep<0) {
      
         MPI_Status  status;
         const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive exahype::records::PatchDescription from node "
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
            msg << "failed to start to receive exahype::records::PatchDescription from node "
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
               msg << "testing for finished receive task for exahype::records::PatchDescription failed: "
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
               "exahype::records::PatchDescription",
               "receive(int)", source,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::PatchDescription",
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
   
   
   
   bool exahype::records::PatchDescription::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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


exahype::records::PatchDescriptionPacked::PersistentRecords::PersistentRecords() {
   
}


exahype::records::PatchDescriptionPacked::PersistentRecords::PersistentRecords(const int& spaceTimePredictorDof, const int& spaceTimeVolumeFluxesDof, const int& dof, const int& updateDof, const int& predictorDof, const int& volumeFluxesDof, const int& extrapolatedPredictorDof, const int& normalFluxesDof, const int& fluctuationsDof, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size):
_spaceTimePredictorDof(spaceTimePredictorDof),
_spaceTimeVolumeFluxesDof(spaceTimeVolumeFluxesDof),
_dof(dof),
_updateDof(updateDof),
_predictorDof(predictorDof),
_volumeFluxesDof(volumeFluxesDof),
_extrapolatedPredictorDof(extrapolatedPredictorDof),
_normalFluxesDof(normalFluxesDof),
_fluctuationsDof(fluctuationsDof),
_level(level),
_offset(offset),
_size(size) {
   
}

exahype::records::PatchDescriptionPacked::PatchDescriptionPacked() {
   
}


exahype::records::PatchDescriptionPacked::PatchDescriptionPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._spaceTimePredictorDof, persistentRecords._spaceTimeVolumeFluxesDof, persistentRecords._dof, persistentRecords._updateDof, persistentRecords._predictorDof, persistentRecords._volumeFluxesDof, persistentRecords._extrapolatedPredictorDof, persistentRecords._normalFluxesDof, persistentRecords._fluctuationsDof, persistentRecords._level, persistentRecords._offset, persistentRecords._size) {
   
}


exahype::records::PatchDescriptionPacked::PatchDescriptionPacked(const int& spaceTimePredictorDof, const int& spaceTimeVolumeFluxesDof, const int& dof, const int& updateDof, const int& predictorDof, const int& volumeFluxesDof, const int& extrapolatedPredictorDof, const int& normalFluxesDof, const int& fluctuationsDof, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size):
_persistentRecords(spaceTimePredictorDof, spaceTimeVolumeFluxesDof, dof, updateDof, predictorDof, volumeFluxesDof, extrapolatedPredictorDof, normalFluxesDof, fluctuationsDof, level, offset, size) {
   
}


exahype::records::PatchDescriptionPacked::~PatchDescriptionPacked() { }



std::string exahype::records::PatchDescriptionPacked::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void exahype::records::PatchDescriptionPacked::toString (std::ostream& out) const {
   out << "("; 
   out << "spaceTimePredictorDof:" << getSpaceTimePredictorDof();
   out << ",";
   out << "spaceTimeVolumeFluxesDof:" << getSpaceTimeVolumeFluxesDof();
   out << ",";
   out << "dof:" << getDof();
   out << ",";
   out << "updateDof:" << getUpdateDof();
   out << ",";
   out << "predictorDof:" << getPredictorDof();
   out << ",";
   out << "volumeFluxesDof:" << getVolumeFluxesDof();
   out << ",";
   out << "extrapolatedPredictorDof:" << getExtrapolatedPredictorDof();
   out << ",";
   out << "normalFluxesDof:" << getNormalFluxesDof();
   out << ",";
   out << "fluctuationsDof:" << getFluctuationsDof();
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


exahype::records::PatchDescriptionPacked::PersistentRecords exahype::records::PatchDescriptionPacked::getPersistentRecords() const {
   return _persistentRecords;
}

exahype::records::PatchDescription exahype::records::PatchDescriptionPacked::convert() const{
   return PatchDescription(
      getSpaceTimePredictorDof(),
      getSpaceTimeVolumeFluxesDof(),
      getDof(),
      getUpdateDof(),
      getPredictorDof(),
      getVolumeFluxesDof(),
      getExtrapolatedPredictorDof(),
      getNormalFluxesDof(),
      getFluctuationsDof(),
      getLevel(),
      getOffset(),
      getSize()
   );
}

#ifdef Parallel
   tarch::logging::Log exahype::records::PatchDescriptionPacked::_log( "exahype::records::PatchDescriptionPacked" );
   
   MPI_Datatype exahype::records::PatchDescriptionPacked::Datatype = 0;
   MPI_Datatype exahype::records::PatchDescriptionPacked::FullDatatype = 0;
   
   
   void exahype::records::PatchDescriptionPacked::initDatatype() {
      {
         PatchDescriptionPacked dummyPatchDescriptionPacked[2];
         
         const int Attributes = 13;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //spaceTimePredictorDof
            MPI_INT,		 //spaceTimeVolumeFluxesDof
            MPI_INT,		 //dof
            MPI_INT,		 //updateDof
            MPI_INT,		 //predictorDof
            MPI_INT,		 //volumeFluxesDof
            MPI_INT,		 //extrapolatedPredictorDof
            MPI_INT,		 //normalFluxesDof
            MPI_INT,		 //fluctuationsDof
            MPI_INT,		 //level
            MPI_DOUBLE,		 //offset
            MPI_DOUBLE,		 //size
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //spaceTimePredictorDof
            1,		 //spaceTimeVolumeFluxesDof
            1,		 //dof
            1,		 //updateDof
            1,		 //predictorDof
            1,		 //volumeFluxesDof
            1,		 //extrapolatedPredictorDof
            1,		 //normalFluxesDof
            1,		 //fluctuationsDof
            1,		 //level
            DIMENSIONS,		 //offset
            DIMENSIONS,		 //size
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._spaceTimePredictorDof))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._spaceTimeVolumeFluxesDof))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._dof))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._updateDof))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._predictorDof))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._volumeFluxesDof))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._extrapolatedPredictorDof))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._normalFluxesDof))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._fluctuationsDof))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._level))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[11] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[1]._persistentRecords._spaceTimePredictorDof))), 		&disp[12] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &PatchDescriptionPacked::Datatype );
         MPI_Type_commit( &PatchDescriptionPacked::Datatype );
         
      }
      {
         PatchDescriptionPacked dummyPatchDescriptionPacked[2];
         
         const int Attributes = 13;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //spaceTimePredictorDof
            MPI_INT,		 //spaceTimeVolumeFluxesDof
            MPI_INT,		 //dof
            MPI_INT,		 //updateDof
            MPI_INT,		 //predictorDof
            MPI_INT,		 //volumeFluxesDof
            MPI_INT,		 //extrapolatedPredictorDof
            MPI_INT,		 //normalFluxesDof
            MPI_INT,		 //fluctuationsDof
            MPI_INT,		 //level
            MPI_DOUBLE,		 //offset
            MPI_DOUBLE,		 //size
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //spaceTimePredictorDof
            1,		 //spaceTimeVolumeFluxesDof
            1,		 //dof
            1,		 //updateDof
            1,		 //predictorDof
            1,		 //volumeFluxesDof
            1,		 //extrapolatedPredictorDof
            1,		 //normalFluxesDof
            1,		 //fluctuationsDof
            1,		 //level
            DIMENSIONS,		 //offset
            DIMENSIONS,		 //size
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._spaceTimePredictorDof))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._spaceTimeVolumeFluxesDof))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._dof))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._updateDof))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._predictorDof))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._volumeFluxesDof))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._extrapolatedPredictorDof))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._normalFluxesDof))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._fluctuationsDof))), 		&disp[8] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._level))), 		&disp[9] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[10] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[11] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyPatchDescriptionPacked[1]._persistentRecords._spaceTimePredictorDof))), 		&disp[12] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &PatchDescriptionPacked::FullDatatype );
         MPI_Type_commit( &PatchDescriptionPacked::FullDatatype );
         
      }
      
   }
   
   
   void exahype::records::PatchDescriptionPacked::shutdownDatatype() {
      MPI_Type_free( &PatchDescriptionPacked::Datatype );
      MPI_Type_free( &PatchDescriptionPacked::FullDatatype );
      
   }
   
   void exahype::records::PatchDescriptionPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      if (communicateSleep<0) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message exahype::records::PatchDescriptionPacked "
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
            msg << "was not able to send message exahype::records::PatchDescriptionPacked "
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
               msg << "testing for finished send task for exahype::records::PatchDescriptionPacked "
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
               "exahype::records::PatchDescriptionPacked",
               "send(int)", destination,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::PatchDescriptionPacked",
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
   
   
   
   void exahype::records::PatchDescriptionPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      if (communicateSleep<0) {
      
         MPI_Status  status;
         const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive exahype::records::PatchDescriptionPacked from node "
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
            msg << "failed to start to receive exahype::records::PatchDescriptionPacked from node "
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
               msg << "testing for finished receive task for exahype::records::PatchDescriptionPacked failed: "
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
               "exahype::records::PatchDescriptionPacked",
               "receive(int)", source,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::PatchDescriptionPacked",
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
   
   
   
   bool exahype::records::PatchDescriptionPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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



