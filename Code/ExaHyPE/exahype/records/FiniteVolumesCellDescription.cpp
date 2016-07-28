#include "exahype/records/FiniteVolumesCellDescription.h"

exahype::records::FiniteVolumesCellDescription::PersistentRecords::PersistentRecords() {
   
}


exahype::records::FiniteVolumesCellDescription::PersistentRecords::PersistentRecords(const int& solverNumber, const double& timeStepSize, const double& timeStamp, const int& solution, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const Type& type):
_solverNumber(solverNumber),
_timeStepSize(timeStepSize),
_timeStamp(timeStamp),
_solution(solution),
_level(level),
_offset(offset),
_size(size),
_type(type) {
   
}

exahype::records::FiniteVolumesCellDescription::FiniteVolumesCellDescription() {
   
}


exahype::records::FiniteVolumesCellDescription::FiniteVolumesCellDescription(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._solverNumber, persistentRecords._timeStepSize, persistentRecords._timeStamp, persistentRecords._solution, persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._type) {
   
}


exahype::records::FiniteVolumesCellDescription::FiniteVolumesCellDescription(const int& solverNumber, const double& timeStepSize, const double& timeStamp, const int& solution, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const Type& type):
_persistentRecords(solverNumber, timeStepSize, timeStamp, solution, level, offset, size, type) {
   
}


exahype::records::FiniteVolumesCellDescription::~FiniteVolumesCellDescription() { }

std::string exahype::records::FiniteVolumesCellDescription::toString(const Type& param) {
   switch (param) {
      case Cell: return "Cell";
   }
   return "undefined";
}

std::string exahype::records::FiniteVolumesCellDescription::getTypeMapping() {
   return "Type(Cell=0)";
}


std::string exahype::records::FiniteVolumesCellDescription::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void exahype::records::FiniteVolumesCellDescription::toString (std::ostream& out) const {
   out << "("; 
   out << "solverNumber:" << getSolverNumber();
   out << ",";
   out << "timeStepSize:" << getTimeStepSize();
   out << ",";
   out << "timeStamp:" << getTimeStamp();
   out << ",";
   out << "solution:" << getSolution();
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
   out << ",";
   out << "type:" << toString(getType());
   out <<  ")";
}


exahype::records::FiniteVolumesCellDescription::PersistentRecords exahype::records::FiniteVolumesCellDescription::getPersistentRecords() const {
   return _persistentRecords;
}

exahype::records::FiniteVolumesCellDescriptionPacked exahype::records::FiniteVolumesCellDescription::convert() const{
   return FiniteVolumesCellDescriptionPacked(
      getSolverNumber(),
      getTimeStepSize(),
      getTimeStamp(),
      getSolution(),
      getLevel(),
      getOffset(),
      getSize(),
      getType()
   );
}

#ifdef Parallel
   tarch::logging::Log exahype::records::FiniteVolumesCellDescription::_log( "exahype::records::FiniteVolumesCellDescription" );
   
   MPI_Datatype exahype::records::FiniteVolumesCellDescription::Datatype = 0;
   MPI_Datatype exahype::records::FiniteVolumesCellDescription::FullDatatype = 0;
   
   
   void exahype::records::FiniteVolumesCellDescription::initDatatype() {
      {
         FiniteVolumesCellDescription dummyFiniteVolumesCellDescription[2];
         
         const int Attributes = 9;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //solverNumber
            MPI_DOUBLE,		 //timeStepSize
            MPI_DOUBLE,		 //timeStamp
            MPI_INT,		 //solution
            MPI_INT,		 //level
            MPI_DOUBLE,		 //offset
            MPI_DOUBLE,		 //size
            MPI_INT,		 //type
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //solverNumber
            1,		 //timeStepSize
            1,		 //timeStamp
            1,		 //solution
            1,		 //level
            DIMENSIONS,		 //offset
            DIMENSIONS,		 //size
            1,		 //type
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStepSize))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStamp))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solution))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._level))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._offset[0]))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._size[0]))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._type))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[1]._persistentRecords._solverNumber))), 		&disp[8] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &FiniteVolumesCellDescription::Datatype );
         MPI_Type_commit( &FiniteVolumesCellDescription::Datatype );
         
      }
      {
         FiniteVolumesCellDescription dummyFiniteVolumesCellDescription[2];
         
         const int Attributes = 9;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //solverNumber
            MPI_DOUBLE,		 //timeStepSize
            MPI_DOUBLE,		 //timeStamp
            MPI_INT,		 //solution
            MPI_INT,		 //level
            MPI_DOUBLE,		 //offset
            MPI_DOUBLE,		 //size
            MPI_INT,		 //type
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //solverNumber
            1,		 //timeStepSize
            1,		 //timeStamp
            1,		 //solution
            1,		 //level
            DIMENSIONS,		 //offset
            DIMENSIONS,		 //size
            1,		 //type
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solverNumber))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStepSize))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._timeStamp))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._solution))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._level))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._offset[0]))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._size[0]))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[0]._persistentRecords._type))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescription[1]._persistentRecords._solverNumber))), 		&disp[8] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &FiniteVolumesCellDescription::FullDatatype );
         MPI_Type_commit( &FiniteVolumesCellDescription::FullDatatype );
         
      }
      
   }
   
   
   void exahype::records::FiniteVolumesCellDescription::shutdownDatatype() {
      MPI_Type_free( &FiniteVolumesCellDescription::Datatype );
      MPI_Type_free( &FiniteVolumesCellDescription::FullDatatype );
      
   }
   
   void exahype::records::FiniteVolumesCellDescription::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      if (communicateSleep<0) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message exahype::records::FiniteVolumesCellDescription "
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
            msg << "was not able to send message exahype::records::FiniteVolumesCellDescription "
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
               msg << "testing for finished send task for exahype::records::FiniteVolumesCellDescription "
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
               "exahype::records::FiniteVolumesCellDescription",
               "send(int)", destination,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::FiniteVolumesCellDescription",
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
   
   
   
   void exahype::records::FiniteVolumesCellDescription::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      if (communicateSleep<0) {
      
         MPI_Status  status;
         const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive exahype::records::FiniteVolumesCellDescription from node "
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
            msg << "failed to start to receive exahype::records::FiniteVolumesCellDescription from node "
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
               msg << "testing for finished receive task for exahype::records::FiniteVolumesCellDescription failed: "
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
               "exahype::records::FiniteVolumesCellDescription",
               "receive(int)", source,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::FiniteVolumesCellDescription",
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
   
   
   
   bool exahype::records::FiniteVolumesCellDescription::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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


exahype::records::FiniteVolumesCellDescriptionPacked::PersistentRecords::PersistentRecords() {
   
}


exahype::records::FiniteVolumesCellDescriptionPacked::PersistentRecords::PersistentRecords(const int& solverNumber, const double& timeStepSize, const double& timeStamp, const int& solution, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const Type& type):
_solverNumber(solverNumber),
_timeStepSize(timeStepSize),
_timeStamp(timeStamp),
_solution(solution),
_level(level),
_offset(offset),
_size(size),
_type(type) {
   
}

exahype::records::FiniteVolumesCellDescriptionPacked::FiniteVolumesCellDescriptionPacked() {
   
}


exahype::records::FiniteVolumesCellDescriptionPacked::FiniteVolumesCellDescriptionPacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._solverNumber, persistentRecords._timeStepSize, persistentRecords._timeStamp, persistentRecords._solution, persistentRecords._level, persistentRecords._offset, persistentRecords._size, persistentRecords._type) {
   
}


exahype::records::FiniteVolumesCellDescriptionPacked::FiniteVolumesCellDescriptionPacked(const int& solverNumber, const double& timeStepSize, const double& timeStamp, const int& solution, const int& level, const tarch::la::Vector<DIMENSIONS,double>& offset, const tarch::la::Vector<DIMENSIONS,double>& size, const Type& type):
_persistentRecords(solverNumber, timeStepSize, timeStamp, solution, level, offset, size, type) {
   
}


exahype::records::FiniteVolumesCellDescriptionPacked::~FiniteVolumesCellDescriptionPacked() { }

std::string exahype::records::FiniteVolumesCellDescriptionPacked::toString(const Type& param) {
   return exahype::records::FiniteVolumesCellDescription::toString(param);
}

std::string exahype::records::FiniteVolumesCellDescriptionPacked::getTypeMapping() {
   return exahype::records::FiniteVolumesCellDescription::getTypeMapping();
}



std::string exahype::records::FiniteVolumesCellDescriptionPacked::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void exahype::records::FiniteVolumesCellDescriptionPacked::toString (std::ostream& out) const {
   out << "("; 
   out << "solverNumber:" << getSolverNumber();
   out << ",";
   out << "timeStepSize:" << getTimeStepSize();
   out << ",";
   out << "timeStamp:" << getTimeStamp();
   out << ",";
   out << "solution:" << getSolution();
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
   out << ",";
   out << "type:" << toString(getType());
   out <<  ")";
}


exahype::records::FiniteVolumesCellDescriptionPacked::PersistentRecords exahype::records::FiniteVolumesCellDescriptionPacked::getPersistentRecords() const {
   return _persistentRecords;
}

exahype::records::FiniteVolumesCellDescription exahype::records::FiniteVolumesCellDescriptionPacked::convert() const{
   return FiniteVolumesCellDescription(
      getSolverNumber(),
      getTimeStepSize(),
      getTimeStamp(),
      getSolution(),
      getLevel(),
      getOffset(),
      getSize(),
      getType()
   );
}

#ifdef Parallel
   tarch::logging::Log exahype::records::FiniteVolumesCellDescriptionPacked::_log( "exahype::records::FiniteVolumesCellDescriptionPacked" );
   
   MPI_Datatype exahype::records::FiniteVolumesCellDescriptionPacked::Datatype = 0;
   MPI_Datatype exahype::records::FiniteVolumesCellDescriptionPacked::FullDatatype = 0;
   
   
   void exahype::records::FiniteVolumesCellDescriptionPacked::initDatatype() {
      {
         FiniteVolumesCellDescriptionPacked dummyFiniteVolumesCellDescriptionPacked[2];
         
         const int Attributes = 9;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //solverNumber
            MPI_DOUBLE,		 //timeStepSize
            MPI_DOUBLE,		 //timeStamp
            MPI_INT,		 //solution
            MPI_INT,		 //level
            MPI_DOUBLE,		 //offset
            MPI_DOUBLE,		 //size
            MPI_INT,		 //type
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //solverNumber
            1,		 //timeStepSize
            1,		 //timeStamp
            1,		 //solution
            1,		 //level
            DIMENSIONS,		 //offset
            DIMENSIONS,		 //size
            1,		 //type
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStepSize))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStamp))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._type))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[1]._persistentRecords._solverNumber))), 		&disp[8] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &FiniteVolumesCellDescriptionPacked::Datatype );
         MPI_Type_commit( &FiniteVolumesCellDescriptionPacked::Datatype );
         
      }
      {
         FiniteVolumesCellDescriptionPacked dummyFiniteVolumesCellDescriptionPacked[2];
         
         const int Attributes = 9;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //solverNumber
            MPI_DOUBLE,		 //timeStepSize
            MPI_DOUBLE,		 //timeStamp
            MPI_INT,		 //solution
            MPI_INT,		 //level
            MPI_DOUBLE,		 //offset
            MPI_DOUBLE,		 //size
            MPI_INT,		 //type
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //solverNumber
            1,		 //timeStepSize
            1,		 //timeStamp
            1,		 //solution
            1,		 //level
            DIMENSIONS,		 //offset
            DIMENSIONS,		 //size
            1,		 //type
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solverNumber))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStepSize))), 		&disp[1] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._timeStamp))), 		&disp[2] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._solution))), 		&disp[3] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._level))), 		&disp[4] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._offset[0]))), 		&disp[5] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._size[0]))), 		&disp[6] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[0]._persistentRecords._type))), 		&disp[7] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyFiniteVolumesCellDescriptionPacked[1]._persistentRecords._solverNumber))), 		&disp[8] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &FiniteVolumesCellDescriptionPacked::FullDatatype );
         MPI_Type_commit( &FiniteVolumesCellDescriptionPacked::FullDatatype );
         
      }
      
   }
   
   
   void exahype::records::FiniteVolumesCellDescriptionPacked::shutdownDatatype() {
      MPI_Type_free( &FiniteVolumesCellDescriptionPacked::Datatype );
      MPI_Type_free( &FiniteVolumesCellDescriptionPacked::FullDatatype );
      
   }
   
   void exahype::records::FiniteVolumesCellDescriptionPacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      if (communicateSleep<0) {
      
         const int result = MPI_Send(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, destination, tag, tarch::parallel::Node::getInstance().getCommunicator());
         if  (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "was not able to send message exahype::records::FiniteVolumesCellDescriptionPacked "
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
            msg << "was not able to send message exahype::records::FiniteVolumesCellDescriptionPacked "
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
               msg << "testing for finished send task for exahype::records::FiniteVolumesCellDescriptionPacked "
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
               "exahype::records::FiniteVolumesCellDescriptionPacked",
               "send(int)", destination,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::FiniteVolumesCellDescriptionPacked",
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
   
   
   
   void exahype::records::FiniteVolumesCellDescriptionPacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise, int communicateSleep) {
      if (communicateSleep<0) {
      
         MPI_Status  status;
         const int   result = MPI_Recv(this, 1, exchangeOnlyAttributesMarkedWithParallelise ? Datatype : FullDatatype, source, tag, tarch::parallel::Node::getInstance().getCommunicator(), &status);
         if ( result != MPI_SUCCESS ) {
            std::ostringstream msg;
            msg << "failed to start to receive exahype::records::FiniteVolumesCellDescriptionPacked from node "
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
            msg << "failed to start to receive exahype::records::FiniteVolumesCellDescriptionPacked from node "
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
               msg << "testing for finished receive task for exahype::records::FiniteVolumesCellDescriptionPacked failed: "
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
               "exahype::records::FiniteVolumesCellDescriptionPacked",
               "receive(int)", source,tag,1
               );
               triggeredTimeoutWarning = true;
            }
            if (
               tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
               (clock()>timeOutShutdown)
            ) {
               tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
               "exahype::records::FiniteVolumesCellDescriptionPacked",
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
   
   
   
   bool exahype::records::FiniteVolumesCellDescriptionPacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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



