package eu.exahype;


import java.io.IOException;

import eu.exahype.analysis.DepthFirstAdapter;
import eu.exahype.node.AAderdgSolver;
import eu.exahype.node.AProject;
import eu.exahype.node.ATwoDimensionalComputationalDomain;
import eu.exahype.node.AThreeDimensionalComputationalDomain;


public class CreateSolverClasses extends DepthFirstAdapter {
    public Boolean valid = true;

    private DirectoryAndPathChecker   _directoryAndPathChecker;

    private String                    _projectName;
    
    private String                    _microarchitecture;

    private java.util.List<String>    _supportedMicroarchitectures;

    private int                       _dimensions;


    public CreateSolverClasses(DirectoryAndPathChecker  directoryAndPathChecker) {
        _directoryAndPathChecker = directoryAndPathChecker;
        _supportedMicroarchitectures = java.util.Arrays.asList("wsm", "snb", "hsw", "knc", "knl", "noarch");
    }

    @Override
    public void inAProject(AProject node) {
        _projectName = node.getName().toString().trim();

        if (node.getSolver().size()==0) { 
            System.out.println( "there are no solvers in the specification file ... nothing to be done" );      
        }
        
        _microarchitecture = node.getArchitecture().toString().trim().toLowerCase();
        if(!_supportedMicroarchitectures.contains(_microarchitecture)) {
            System.out.println( "Unknown architecture specified ... fallback solution \"noarch\" taken" );
            _microarchitecture = "noarch";
        }

    } 


    @Override
    public void inATwoDimensionalComputationalDomain(ATwoDimensionalComputationalDomain node) {
        _dimensions = 2;
    }


    @Override
    public void inAThreeDimensionalComputationalDomain(AThreeDimensionalComputationalDomain node) {
        _dimensions = 3;
    }


    @Override
    public void inAAderdgSolver(AAderdgSolver node) {
        String solverName = node.getName().toString().trim();

        java.io.File headerFile                  = new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + ".h");
        java.io.File userImplementationFile      = new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + ".cpp");
        java.io.File generatedImplementationFile = new java.io.File(_directoryAndPathChecker.outputDirectory.getAbsolutePath() + "/" + solverName + "_generated.cpp");

        String kernel = node.getKernel().toString().trim();

        try {      
            // =====================
            // Write all the headers
            // =====================
            if (headerFile.exists()) {
                System.out.println( "create header of solver " + solverName + " ... header " + headerFile.getAbsoluteFile() + " does exist already. Remove to allow toolkit to regenerate it (changes will be lost)" );      
            }
            else {
                java.io.BufferedWriter headerWriter = new java.io.BufferedWriter(new java.io.FileWriter(headerFile));

                if (kernel.equals("user::defined")) {
                    System.out.println( "create header of solver " + solverName + " ... ok" );      
                    writeMinimalADERDGSolverHeader( solverName, headerWriter );
                }
                else if (kernel.equals("user::fluxes")) {
                    System.out.println( "create header of solver " + solverName + " ... ok" );      

                    writeADERDGSolverHeaderForUserFluxes( solverName, headerWriter );
                }
                else if (node.getKernel().toString().trim().equals("kernel::euler2d")) {
                    System.out.println( "create header of solver " + solverName + " ... ok" );      
                    writeMinimalADERDGSolverHeader( solverName, headerWriter );
                }
                else {
                    System.err.println( "ERROR: unknown ADER-DG kernel type " + kernel + " ... failed" );      
                    valid = false;
                }

                headerWriter.close();
            }

            // =======================================
            // Write all the user implementation files
            // =======================================
            if (userImplementationFile.exists()) {
                System.out.println( "user's implementation file of solver " + solverName + " ... does exist already. Is not overwritten" );      
            }
            else {
                java.io.BufferedWriter userImplementationWriter = new java.io.BufferedWriter(new java.io.FileWriter(userImplementationFile));

                if (kernel.equals("user::defined")) {
                    System.out.println( "create user implementation template of solver " + solverName + " ... please complete" );      

                    writeADERDGSolverUserImplementationForUserDefined(solverName, node.getVariables().toString().trim(), node.getOrder().toString().trim(), userImplementationWriter);
                }
                else if (kernel.equals("user::fluxes")) {
                    System.out.println( "create user implementation template of solver " + solverName + " ... please complete" );      

                    writeADERDGSolverUserImplementationForUserFluxes(solverName, node.getVariables().toString().trim(), node.getOrder().toString().trim(), userImplementationWriter);
                }
                else if (node.getKernel().toString().trim().equals("kernel::euler2d")) {
                    System.out.println( "create user implementation template of solver " + solverName + " ... please complete" );      

                    writeMinimalADERDGSolverUserImplementation(solverName, node.getVariables().toString().trim(), node.getOrder().toString().trim(), userImplementationWriter);
                }
                else {
                    System.err.println( "ERROR: unknown ADER-DG kernel type " + kernel + " ... failed" );      
                    valid = false;
                }

                userImplementationWriter.close();
            }


            // ============================================
            // Write all the generated implementation files
            // ============================================
            if (generatedImplementationFile.exists()) {
                System.out.println( "generated implementation file of solver " + solverName + " ... does exist already. Is overwritten" );      
            }
            java.io.BufferedWriter generatedImplementationWriter = new java.io.BufferedWriter(new java.io.FileWriter(generatedImplementationFile));

            if (kernel.equals("user::defined")) {
                System.out.println( "create generated implementation of solver " + solverName + " ... is empty (ok)" );      

                writeADERDGSolverGeneratedImplementationForUserDefined(solverName, generatedImplementationWriter);
            }
            else if (kernel.equals("user::fluxes")) {
                System.out.println( "create generated implementation of solver " + solverName + " ... ok" );      

                writeADERDGSolverGeneratedImplementationForUserFluxes(solverName, node.getVariables().toString().trim(), node.getOrder().toString().trim(), generatedImplementationWriter);
            }
            else if (node.getKernel().toString().trim().equals("kernel::euler2d")) {
                System.out.println( "create generated implementation of solver " + solverName + " ... ok" );     

                writeADERDGSolverGeneratedImplementationForKernelEuler2d(solverName, node.getVariables().toString().trim(), node.getOrder().toString().trim(), generatedImplementationWriter);
                invokeCodeGenerator(solverName, node.getVariables().toString().trim(), node.getOrder().toString().trim());
            }
            else {
                System.err.println( "ERROR: unknown ADER-DG kernel type " + kernel + " ... failed" );      
                valid = false;
            }

            generatedImplementationWriter.close();
        }      
        catch (Exception exc) {
            System.err.println( "ERROR: " + exc.toString() );
            valid = false;
        }
    }


    /**
     * Write header with ExaHyPE copyright. Should be inserted for any solver's 
     * header.
     */
    private void writeHeaderCopyright( java.io.BufferedWriter writer ) throws IOException {
        writer.write("// This file is generated by the ExaHyPE toolkit.\n");
        writer.write("// Please do not modify - it will be overwritten by the next\n");
        writer.write("// ExaHyPE toolkit call.\n");
        writer.write("// \n");
        writer.write("// ========================\n");
        writer.write("//   www.exahype.eu\n");
        writer.write("// ========================\n");
    }

    /**
     * Adds all the default includes of any solver as well as the solver define. 
     * Is used by all solvers.
     */
    private void writeHeaderIncludesAndDefines( java.io.BufferedWriter writer, String solverName ) throws IOException {
        writer.write("\n\n\n");
        writer.write("#include \"exahype/solvers/Solver.h\"");
        writer.write("\n\n\n");

        writer.write("namespace " + _projectName + "{\n");
        writer.write("  class " + solverName + ";\n");
        writer.write("}\n\n\n");
    }

    /**
     * Creates all the public operations that are mandatory for any solver.
     */
    private void writeHeaderMinimalADERDGClassSignature( java.io.BufferedWriter writer, String solverName ) throws IOException {
        writer.write("class " + _projectName + "::" + solverName + ": public exahype::solvers::Solver {\n");
        writer.write("  public:\n");
        writer.write("    " + solverName + "(int kernelNumber); \n");
        writer.write("    virtual int getMinimumTreeDepth() const;\n");

        writer.write("    virtual void spaceTimePredictor(double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ); \n");
        writer.write("    virtual void solutionUpdate(double* luh, const double* const lduh, const double dt);\n");
        writer.write("    virtual void volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx);\n" );
        writer.write("    virtual void surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx);\n" );
        writer.write("    virtual void riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex);\n" );
        writer.write("    virtual double stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx );\n" );
        writer.write("    virtual void initialCondition(double* luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx);\n" );
    }

    private void writeMinimalADERDGSolverHeader(
            String                 solverName,
            java.io.BufferedWriter writer
            ) throws IOException {
        writeHeaderCopyright(writer);
        writeHeaderIncludesAndDefines(writer,solverName);
        writeHeaderMinimalADERDGClassSignature(writer,solverName);

        writer.write("};\n\n\n");
    }


    private void writeADERDGSolverHeaderForUserFluxes(
            String                 solverName,
            java.io.BufferedWriter writer
            ) throws IOException {
        writeHeaderCopyright(writer);
        writeHeaderIncludesAndDefines(writer,solverName);
        writeHeaderMinimalADERDGClassSignature(writer,solverName);

        writer.write("  private:\n");
        if (_dimensions==2) {
          writer.write("    static void flux(const double* const Q, double* f, double* g);\n");
        }
        else {
         writer.write("    static void flux(const double* const Q, double* f, double* g, double* h);\n");
        }
        writer.write("    static void eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda);\n");
        writer.write("    static void initialValues(const double* const x, double* Q);\n" );

        writer.write("};\n\n\n");
    }


    private void writeADERDGSolverGeneratedImplementationForUserDefined(
            String                 solverName,
            java.io.BufferedWriter writer
            ) throws IOException {
        writer.write("// ==============================================\n");
        writer.write("// Please do not change the implementations below\n");
        writer.write("// =============================---==============\n");
        writer.write("#include \"" + solverName + ".h\"\n");
        writer.write("\n\n\n");
        writer.write("// This file is empty as a user::defined kernel is chosen, i.e. the user\n");
        writer.write("// wants to implement everything.");
        writer.write("\n\n\n");
    }


    private void writeADERDGSolverGeneratedImplementationForUserFluxes(
            String                 solverName,
            String                 numberOfVariables,
            String                 order,
            java.io.BufferedWriter writer
            ) throws IOException {
        writer.write("// ==============================================\n");
        writer.write("// Please do not change the implementations below\n");
        writer.write("// =============================---==============\n");
        writer.write("#include \"" + solverName + ".h\"\n");
        writer.write("#include \"kernels/aderdg/generic/Kernels.h\"\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::spaceTimePredictor( double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ) {\n");
        writer.write("   kernels::aderdg::generic::spaceTimePredictor<flux>( lQi, lFi, lQhi, lFhi, lQhbnd, lFhbnd, luh, dx, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::solutionUpdate(double* luh, const double* const lduh, const double dt) {\n");
        writer.write("   kernels::aderdg::generic::solutionUpdate( luh, lduh, dt, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
        writer.write("   kernels::aderdg::generic::volumeIntegral( lduh, lFhi, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
        writer.write("   kernels::aderdg::generic::surfaceIntegral( lduh, lFhbnd, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) {\n");
        writer.write("   kernels::aderdg::generic::riemannSolver<eigenvalues>( FL, FR, QL, QR, dt, normalNonZeroIndex, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "double " + _projectName + "::" + solverName + "::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
        writer.write("   return kernels::aderdg::generic::stableTimeStepSize<eigenvalues>( luh, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::initialCondition(double* luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
        writer.write("   kernels::aderdg::generic::initialCondition<initialValues>( luh, center, dx, getNumberOfVariables(), getNodesPerCoordinateAxis() );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
    }


    private void writeADERDGSolverGeneratedImplementationForKernelEuler2d(
            String                 solverName,
            String                 numberOfVariables,
            String                 order,
            java.io.BufferedWriter writer
            ) throws IOException {
        writer.write("// ==============================================\n");
        writer.write("// Please do not change the implementations below\n");
        writer.write("// =============================---==============\n");
        writer.write("#include \"" + solverName + ".h\"\n");
        writer.write( "#include \"kernels/aderdg/optimised/defines.h\"\n");
        writer.write( "#include \"kernels/aderdg/optimised/Kernels.h\"\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::spaceTimePredictor( double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ) {\n");
        writer.write("   kernels::aderdg::optimised::picardLoop<flux>( lQi, lFi, luh, dx, dt );\n");
        writer.write("   kernels::aderdg::optimised::predictor( lQhi, lFhi, lQi, lFi );\n");
        writer.write("   kernels::aderdg::optimised::extrapolator( lQhbnd, lFhbnd, lQhi, lFhi );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::solutionUpdate(double* luh, const double* const lduh, const double dt) {\n");
        writer.write("   kernels::aderdg::optimised::solutionUpdate( luh, lduh, dt );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
        writer.write("   kernels::aderdg::optimised::volumeIntegral( lduh, lFhi, dx );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
        writer.write("   kernels::aderdg::optimised::surfaceIntegral( lduh, lFhbnd, dx );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) {\n");
        writer.write("   kernels::aderdg::optimised::riemannSolver<eigenvalues>( FL, FR, QL, QR, dt, normalNonZeroIndex );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "double " + _projectName + "::" + solverName + "::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
        writer.write("   return kernels::aderdg::optimised::stableTimeStepSize<eigenvalues>( luh, dx );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write( "void " + _projectName + "::" + solverName + "initialCondition(double* luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
        writer.write("   kernels::aderdg::optimised::initialCondition<initialValues>( luh, center, dx );\n");
        writer.write("}\n");
        writer.write("\n\n\n");
    }


    private void writeADERDGSolverUserImplementationForUserFluxes(
            String                 solverName,
            String                 numberOfVariables,
            String                 order,
            java.io.BufferedWriter writer
            ) throws IOException {
        writeMinimalADERDGSolverUserImplementation(solverName,numberOfVariables,order,writer);

            int digits = String.valueOf(numberOfVariables).length();

            if (_dimensions==2) {
              writer.write("void " + _projectName + "::" + solverName + "::flux(const double* const Q, double* f, double* g) {\n");
            }
            else {
              writer.write("void " + _projectName + "::" + solverName + "::flux(const double* const Q, double* f, double* g, double* h) {\n");
            }
            writer.write("  // Dimensions             = "+_dimensions      +"\n");
            writer.write("  // Number of variables    = "+numberOfVariables+"\n");
            writer.write("  // f\n");
            writer.write("  // @todo Please implement\n");
            for (int i=0; i < Integer.parseInt(numberOfVariables); i++) {
                writer.write("  f["+String.format("%"+digits+"d",i)+"] = 0.0;\n");    
            }
            writer.write("  // g\n");
            writer.write("  // @todo Please implement\n");
            for (int i=0; i < Integer.parseInt(numberOfVariables); i++) {
                writer.write("  g["+String.format("%"+digits+"d",i)+"] = 0.0;\n");    
            }
            if (_dimensions==3) {
              writer.write("  // h\n");
              writer.write("  // @todo Please implement\n");
              for (int i=0; i < Integer.parseInt(numberOfVariables); i++) {
                  writer.write("  h["+String.format("%"+digits+"d",i)+"] = 0.0;\n");    
              }
            }
            writer.write("}\n");
            writer.write("\n\n\n");
            writer.write("void " + _projectName + "::" + solverName + "::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {\n");
            writer.write("  // Dimensions             = "+_dimensions      +"\n");
            writer.write("  // Number of variables    = "+numberOfVariables+"\n");
            writer.write("  // @todo Please implement\n");
            for (int i=0; i < Integer.parseInt(numberOfVariables); i++) {
                writer.write("  lambda["+String.format("%"+digits+"d",i)+"] = 0.0;\n");    
            }
            writer.write("}\n");
            writer.write("\n\n\n");
            writer.write("void " + _projectName + "::" + solverName + "::initialValues(const double* const x, double* Q) {\n");
            writer.write("  // Dimensions             = "+_dimensions      +"\n");
            writer.write("  // Number of variables    = "+numberOfVariables+"\n");
            writer.write("  // @todo Please implement\n");
            for (int i=0; i < Integer.parseInt(numberOfVariables); i++) {
                writer.write("  Q["+String.format("%"+digits+"d",i)+"] = 0.0;\n");    
            }
            writer.write("}\n");
            writer.write("\n\n\n");
    }


    private void writeADERDGSolverUserImplementationForUserDefined(
            String                 solverName,
            String                 numberOfVariables,
            String                 order,
            java.io.BufferedWriter writer
            ) throws IOException {
        writer.write("#include \"" + solverName + ".h\"\n");
        writer.write("\n\n\n");
        writer.write(_projectName + "::" + solverName + "::" + solverName + "( int kernelNumber):\n");
        writer.write("  exahype::solvers::Solver(\"" + solverName + "\",exahype::solvers::Solver::ADER_DG,kernelNumber," + numberOfVariables + "," + order + "+1) {\n");
        writer.write("  // @todo Please implement/augment if required\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write("int " + _projectName + "::" + solverName + "::getMinimumTreeDepth() const {\n");
        writer.write("  // @todo Please implement\n");
        writer.write("  return 3;\n");
        writer.write("}\n");
        writer.write("\n\n\n");

            writer.write("void " + _projectName + "::" + solverName + "::spaceTimePredictor(double* lQi, double* lFi, double* lQhi, double* lFhi, double* lQhbnd, double* lFhbnd, const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx, const double dt ) {\n");
            writer.write("  // @todo Please implement\n");
            writer.write("}\n");
            writer.write("\n\n\n");
            writer.write("void " + _projectName + "::" + solverName + "::solutionUpdate(double* luh, const double* const lduh, const double dt) {\n");
            writer.write("  // @todo Please implement\n");
            writer.write("}\n");
            writer.write("\n\n\n");
            writer.write("void " + _projectName + "::" + solverName + "::volumeIntegral(double* lduh, const double* const lFhi, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
            writer.write("  // @todo Please implement\n");
            writer.write("}\n");
            writer.write("\n\n\n");
            writer.write("void " + _projectName + "::" + solverName + "::surfaceIntegral(double* lduh, const double* const lFhbnd, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
            writer.write("  // @todo Please implement\n");
            writer.write("}\n");
            writer.write("\n\n\n");
            writer.write("void " + _projectName + "::" + solverName + "::riemannSolver(double* FL, double* FR, const double* const QL, const double* const QR, const double dt, const int normalNonZeroIndex) {\n");
            writer.write("  // @todo Please implement\n");
            writer.write("}\n");
            writer.write("\n\n\n");
            writer.write("double " + _projectName + "::" + solverName + "::stableTimeStepSize(const double* const luh, const tarch::la::Vector<DIMENSIONS,double>& dx ) {\n");
            writer.write("  // @todo Please implement\n");
            writer.write("  return 1.0;\n");
            writer.write("}\n");
            writer.write("\n\n\n");
            writer.write("void " + _projectName + "::" + solverName + "::initialCondition(double* luh, const tarch::la::Vector<DIMENSIONS,double>& center, const tarch::la::Vector<DIMENSIONS,double>& dx) {\n");
            writer.write("  // @todo Please implement\n");
            writer.write("}\n");
            writer.write("\n\n\n");
    }


    private void writeMinimalADERDGSolverUserImplementation(
            String                 solverName,
            String                 numberOfVariables,
            String                 order,
            java.io.BufferedWriter writer
            ) throws IOException {
        writer.write("#include \"" + solverName + ".h\"\n");
        writer.write("\n\n\n");
        writer.write(_projectName + "::" + solverName + "::" + solverName + "( int kernelNumber):\n");
        writer.write("  exahype::solvers::Solver(\"" + solverName + "\",exahype::solvers::Solver::ADER_DG,kernelNumber," + numberOfVariables + "," + order + "+1) {\n");
        writer.write("  // @todo Please implement/augment if required\n");
        writer.write("}\n");
        writer.write("\n\n\n");
        writer.write("int " + _projectName + "::" + solverName + "::getMinimumTreeDepth() const {\n");
        writer.write("  // @todo Please implement\n");
        writer.write("  return 3;\n");
        writer.write("}\n");
    }

    private void invokeCodeGenerator(
            String                 solverName,
            String                 numberOfVariables,
            String                 order
            ) throws IOException {
        // TODO adapt path
        String currentDirectory = System.getProperty("user.dir");
        java.nio.file.Path pathToCodeGenerator = java.nio.file.Paths.get(currentDirectory+"/Miscellaneous/CodeGenerator/Driver.py");
		    if(java.nio.file.Files.notExists(pathToCodeGenerator)) {
			    System.err.println("ERROR: Code generator not found. Can't generated optimised kernels.");
			    return;
		    }

		    // set up the command to execute the code generator
		    String args          = " " + solverName                    + " "
		                               + numberOfVariables             + " "
		                               + order                         + " "
		                               + Integer.toString(_dimensions) + " "
		                               + _microarchitecture            + " "
		                               + "--precision=DP";  //double precision
		                               
		    String bashCommand   = "python " + pathToCodeGenerator + args ;

		    Runtime runtime = Runtime.getRuntime();

	      // execute the command line program
	      Process codeGenerator = runtime.exec(bashCommand);

	      // capture any output that is produced by the code generator and print it line-by-line
	      java.io.BufferedReader codeGeneratorsOutputReader = new java.io.BufferedReader(new java.io.InputStreamReader(codeGenerator.getInputStream()));
	      String line = "";
	      while((line = codeGeneratorsOutputReader.readLine()) != null) {
		        System.out.println(line);
	      }
    }
} 
