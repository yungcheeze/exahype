package eu.exahype.variables;

import java.io.BufferedWriter;
import java.util.LinkedHashMap;
import java.util.Map;

import eu.exahype.IOUtils;
import eu.exahype.node.AAderdgSolver;
import eu.exahype.node.AFiniteVolumesSolver;
import eu.exahype.node.ALimitingAderdgSolver;
import eu.exahype.node.PSolver;

public class Variables {
  int                 _dimensions;
  Map<String,Integer> _variablesMap;
  int                 _numberOfVariables;
  Map<String,Integer> _parametersMap;
  int                 _numberOfParameters;
  
  public Map<String, Integer> get_variablesMap() {
    return _variablesMap;
  }

  public int get_numberOfVariables() {
    return _numberOfVariables;
  }

  public Map<String, Integer> get_parametersMap() {
    return _parametersMap;
  }

  public int get_numberOfParameters() {
    return _numberOfParameters;
  }

  /**
   * @note We rely on a linked hash map here that does not change the order of the variables.
   */
  public static Map<String, Integer> readVariables(PSolver node) {
    String variablesAsString = null;
    
    if (node instanceof AAderdgSolver) {
      variablesAsString = ((AAderdgSolver) node).getVariables().getText();
    } else if (node instanceof ALimitingAderdgSolver) {
      variablesAsString = ((ALimitingAderdgSolver) node).getVariables().getText();
    } else if (node instanceof AFiniteVolumesSolver) {
      variablesAsString = ((AFiniteVolumesSolver) node).getVariables().getText();
    } else {
      System.out.println("ERROR: I do not know how to handle solver type "+node.getClass().toString()+"!");
      System.exit(1);
    }
    
    try { // the user only gave us a number, e.g., 5, instead of a list, e.g., v0:1, v1:3, v2:3.
      int numberOfVariables = Integer.parseInt(variablesAsString);
      
      Map<String, Integer> map = new LinkedHashMap<String, Integer>(1);
      map.put("Q", numberOfVariables);
      return map;
    } catch (NumberFormatException exception) { // the user gave us a list       
      String[] variables = variablesAsString.split(",");
      Map<String, Integer> map = new LinkedHashMap<String, Integer>(variables.length);
      
      for (String variable : variables) {
        String[] identifierAndQuantity = variable.split(":");
        
        String identifier = identifierAndQuantity[0].trim();
        try {
          int multiplicity    = Integer.parseInt(identifierAndQuantity[1].trim());
          
          if (multiplicity <= 0) {
            System.out.println("ERROR: Quantity specifier of '"+identifier+"' is not a positive integer!");
            System.exit(1);
          }
          
          map.put(identifier, multiplicity);
          
          // System.out.println("Found variable "+identifier+" with "+multiplicity+" elements."); // Comment in for debugging purposes
          
        } catch (NumberFormatException exception2) { 
          System.out.println("ERROR: Quantity specifier of '"+identifier+"' is not a positive integer!");
          System.exit(1);
        }
      }
      return map;
    }
  }
  
  /**
   * @note We rely on a linked hash map here that does not change the order of the parameters.
   */
  public static Map<String, Integer> readParameters(PSolver node) {
    String parametersAsString = null;
    
    if (node instanceof AAderdgSolver) {
      parametersAsString = ((AAderdgSolver) node).getParameters().getText();
    } else if (node instanceof ALimitingAderdgSolver) {
      parametersAsString = ((ALimitingAderdgSolver) node).getParameters().getText();
    } else if (node instanceof AFiniteVolumesSolver) {
      parametersAsString = ((AFiniteVolumesSolver) node).getParameters().getText();
    } else {
      System.out.println("ERROR: I do not know how to handle solver type "+node.getClass().toString()+"!");
      System.exit(1);
    }
    
    try { // the user only gave us a number, e.g., 5, instead of a list, e.g., v0:1, v1:3, v2:3.
      int numberOfVariables = Integer.parseInt(parametersAsString);
      
      Map<String, Integer> map = new LinkedHashMap<String, Integer>(1);
      map.put("Q", numberOfVariables);
      return map;
    } catch (NumberFormatException exception) { // the user gave us a list       
      String[] variables = parametersAsString.split(",");
      Map<String, Integer> map = new LinkedHashMap<String, Integer>(variables.length);
      
      for (String variable : variables) {
        String[] identifierAndQuantity = variable.split(":");
        
        String identifier = identifierAndQuantity[0].trim();
        try {
          int dimension    = Integer.parseInt(identifierAndQuantity[1].trim());
          
          if (dimension <= 0) {
            System.out.println("ERROR: Quantity specifier of '"+identifier+"' is not a positive integer!");
            System.exit(1);
          }
          
          map.put(identifier, dimension);
          
          // System.out.println("Found variable "+identifier+" with "+dimension+" elements."); // Comment in for debugging purposes
          
        } catch (NumberFormatException exception2) { 
          System.out.println("ERROR: Quantity specifier of '"+identifier+"' is not a positive integer!");
          System.exit(1);
        }
      }
      return map;
    }
  }
  
  public static int sumMultiplicities(Map<String,Integer> variables) {
    int sumOfMultiplicities = 0;
    for (String key : variables.keySet()) {
      sumOfMultiplicities += variables.get(key);
    }
    
    return sumOfMultiplicities;
  }
  
  public Variables(PSolver node, int dimensions) {
    _variablesMap       = readVariables(node);
    _parametersMap      = readParameters(node);
    _numberOfVariables  = sumMultiplicities(_variablesMap);
    _numberOfParameters = sumMultiplicities(_parametersMap);
    _dimensions         = dimensions;
  }
  
  public Variables(Map<String, Integer> variablesMap, Map<String, Integer> parametersMap, int dimensions) {
    _variablesMap       = variablesMap;
    _parametersMap      = parametersMap;
    _numberOfVariables  = sumMultiplicities(variablesMap);
    _numberOfParameters = sumMultiplicities(parametersMap);
    _dimensions         = dimensions;
  }
  
  private String appendVariableGetters(String getters, String identifier, int multiplicity, int offset) {
    String indent  = "    ";
    
    assert multiplicity > 0;
    if (multiplicity==1) {
      // ex: double rho() const;
      getters += indent + "double "+identifier+"() const { return _Q["+offset+"]; }\n\n";
    } else {
      // ex: double v(index) const;
      getters += indent + "double "+identifier+"(int index) const {\n"
              +  indent + "  assertion(index >= 0 && index<"+multiplicity+");\n"
              +  indent + "  return _Q["+offset+"+index];\n"
              +  indent + "}\n\n";
      // ex: tarch::la::Vector<3,double> v() const;
      getters += indent + "tarch::la::Vector<"+multiplicity+",double> "+identifier+"() const {\n"
               + indent + "  tarch::la::Vector<"+multiplicity+",double> values(_Q+"+offset+");\n"
               + indent + "  return values;\n"
               + indent + "}\n\n";
    }
    return getters;
  }
  
  private String createVariablesGetters() {
    String getters = "";
    
    int offset = 0;
    for (String identifier : _variablesMap.keySet()) {
      int multiplicity = _variablesMap.get(identifier);
      getters  = appendVariableGetters(getters,identifier,multiplicity,offset);
      offset  += multiplicity;
    }
    for (String identifier : _parametersMap.keySet()) {
      int multiplicity = _parametersMap.get(identifier);
      getters  = appendVariableGetters(getters,identifier,multiplicity,offset);
      offset  += multiplicity;
    }
    
    return getters;
  }
  
  private String appendVariableSetters(String setters, String identifier, int multiplicity, int offset) {
    String indent  = "    ";
    
    assert multiplicity > 0;
    if (multiplicity==1) {
      setters += indent + "double& "+identifier+"() { return _Q["+offset+"]; }\n\n";
    } else {
      // ex: double& v(int index);
      setters += indent +"double& "+identifier+"(int index) { return _Q["+offset+"+index]; }\n\n";
      
      // ex: void v(const tarch::la::Vector<3,double>& values);
      setters += indent +"void "+identifier+"(const tarch::la::Vector<"+multiplicity+",double>& values) {\n";
      for (int i=0; i<multiplicity; i++) {
        setters += indent +"  "+"*(_Q+"+(offset+i)+")=values["+i+"];\n";
      }
      setters += indent +"}\n\n";
      
      // ex: void v(double v0, double v1, double v2);
      setters += indent +"void "+identifier+"(";
      for (int i=0; i<multiplicity; i++) {
        setters += identifier+i+( i<multiplicity-1 ? "," : ") {\n" );
      }
      for (int i=0; i<multiplicity; i++) {
        setters += indent +"  "+"*(_Q+"+(offset+i)+")="+identifier+i+";\n";
      }
      setters += indent + "}\n\n";
    }
    return setters;
  }
  
  private String createVariablesSetters() {
    String setters = "";
    
    int offset = 0;
    for (String identifier : _variablesMap.keySet()) {
      int multiplicity = _variablesMap.get(identifier);
      setters = appendVariableSetters(setters, identifier, multiplicity, offset);
      offset += multiplicity;
    }
    for (String identifier : _parametersMap.keySet()) {
      int multiplicity = _parametersMap.get(identifier);
      setters = appendVariableSetters(setters, identifier, multiplicity, offset);
      offset += multiplicity;
    }
    
    return setters;
  }
  
  /**
   * @note The current implementation assumes a column major memory layout
   * of the fluxes.
   */
  private String appendFluxesGetters(String getters, String identifier,int multiplicity, int offset) {
    String indent  = "    ";
    
    assert multiplicity > 0;
    if (multiplicity==1) {
      // ex: double rho(int column) const;  
      getters += indent + "double "+identifier+"(int column) const {\n"
              +  indent + "  assertion(column >= 0 && column<DIMENSIONS);\n"
              +  indent + "  return _F[column]["+offset+"+row];\n"
              +  indent + "}\n\n";
      // ex: tarch::la::Vector<DIMENSIONS,double> rho() const;
      getters += indent + "tarch::la::Vector<DIMENSIONS,double> "+identifier+"() const {\n"
              +  indent + "  tarch::la::Vector<DIMENSIONS,double> values(";
              for (int i=0; i<_dimensions; i++) {
                getters += "_F["+i+"]["+offset+"]" + ( i<_dimensions-1 ? "," : ");\n" );
              }
      getters += indent + "  return values;\n"
              +  indent + "}\n\n";
      
    } else {
      // ex: double v(int row, int column) const;  
      getters += indent + "double "+identifier+"(int row, int column) const {\n"
              +  indent + "  assertion(column >= 0 && column<DIMENSIONS);\n"
              +  indent + "  assertion(row >= 0 && row<"+multiplicity+");\n"
              +  indent + "  return _F[column]["+offset+"+row];\n"
              +  indent + "}\n\n";
      
      // tarch::la::Vector<3,double> v(int row) const;
      getters += indent + "tarch::la::Vector<DIMENSIONS,double> "+identifier+"(int row) const {\n"
              +  indent + "  tarch::la::Vector<DIMENSIONS,double> values(";
      for (int i=0; i<_dimensions; i++) {
        getters += "_F["+i+"]["+offset+"+row]" + ( i<_dimensions-1 ? "," : ");\n" );
      }
      getters += indent + "  return values;\n"
              +  indent + "}\n\n";
      
      // ex: tarch::la::Matrix<3,3,double> v() const;
      getters += indent + "tarch::la::Matrix<"+multiplicity+",DIMENSIONS,double> "+identifier+"() const {\n"
              +  indent + "  tarch::la::Matrix<"+multiplicity+",DIMENSIONS,double> values;\n";
      getters += indent + "  values = ";      
      for (int i=0; i<multiplicity; i++) {
        for (int j=0; j<_dimensions; j++) {
          getters += "_F["+j+"]["+(offset+i)+"]" + ( j<_dimensions-1 ? "," : "" );
        }
        getters += ( i<multiplicity-1 ? ",\n"+indent + "           " : ";\n" );
      }
      getters += indent + "  return values;\n"
              +  indent + "}\n\n";
      
    }
    return getters;
  }
  
  private String createFluxesGetters() {
    String getters = "";
    
    int offset = 0;
    for (String identifier : _variablesMap.keySet()) {
      int multiplicity = _variablesMap.get(identifier);
      getters  = appendFluxesGetters(getters,identifier,multiplicity,offset);
      offset  += multiplicity;
    }
    
    return getters;
  }

  private String appendFluxesSetters(String setters, String identifier,
      int multiplicity, int offset) {
    String indent  = "    ";

    assert multiplicity > 0;
    if (multiplicity==1) {
      // ex: double& v(int column);
      setters += indent + "double& "+identifier+"(int column) {\n" 
              +  indent + "  assertion(column >= 0 && column<DIMENSIONS);\n"
              +  indent + "  return _F[column]["+offset+"];\n" 
              +  indent + "}\n\n";
      
      // ex: void rho(tarch::la::Vector<DIMENSIONS,double>& values);
      setters += indent +"void "+identifier+"(tarch::la::Vector<DIMENSIONS,double>& values) {\n";
      for (int j=0; j<_dimensions; j++) {
        setters += indent +"  "+"_F["+j+"]["+offset+"]=values["+j+"];\n";
      }
      setters += indent + "}\n\n";
      
      // ex: void v(double v0, double v1, double v2);
      setters += indent +"void "+identifier+"(";
      for (int j=0; j<_dimensions; j++) {
        setters += "v"+j+( j<_dimensions-1 ? "," : ") {\n" );
      }
      for (int j=0; j<_dimensions; j++) {
        setters += indent +"  "+"_F["+j+"]["+offset+"]=v"+j+";\n";
      }
      setters += indent + "}\n\n";
      
    } else {
      // ex: double& v(int row, int column);
      setters += indent + "double& "+identifier+"(int row, int column) {\n"
              +  indent + "  assertion(row >= 0 && row<"+multiplicity+");\n"
              +  indent + "  assertion(column >= 0 && column<DIMENSIONS);\n"
              +  indent + "  return _F[column]["+offset+"+row];\n"
              +  indent + "}\n\n";

      // ex: void v(int row, const tarch::la::Vector<3,double>& values);
      setters += indent +"void "+identifier+"(int row, const tarch::la::Vector<DIMENSIONS,double>& values) {\n";
      for (int j=0; j<_dimensions; j++) {
        setters += indent +"  "+"_F["+j+"]["+offset+"+row]=values["+j+"];\n";
      }
      setters += indent +"}\n\n";
      
      // ex: void v(const tarch::la::Matrix<3,DIMENSIONS,double>& values);
      setters += indent +"void "+identifier+"(int row, const tarch::la::Matrix<"+multiplicity+",DIMENSIONS,double>& values) {\n";
      for (int i=0; i<multiplicity; i++) {
        for (int j=0; j<_dimensions; j++) {
          setters += indent +"  "+"_F["+j+"]["+(offset+i)+"]=values["+i+"]["+j+"];\n";
        }
      }
      setters += indent +"}\n\n";
      
      // ex: void v(int row, double v0, double v1, double v2);
      setters += indent +"void "+identifier+"(int row, ";
      for (int j=0; j<_dimensions; j++) {
        setters += "int v"+j+( j<_dimensions-1 ? "," : ") {\n" );
      }
      for (int j=0; j<_dimensions; j++) {
        setters += indent +"  "+"_F["+j+"]["+offset+"+row]=v"+j+";\n";
      }
      setters += indent +"}\n\n";

      // ex: void v(double v00, double v01, double v02, double v10, double v11, double v12, double v20, double v21, double v22);
      setters += indent +"void "+identifier+"(";
      for (int i=0; i<multiplicity; i++) {
        for (int j=0; j<_dimensions; j++) {
          setters += "int v"+i+j+( j<_dimensions-1 ? "," : "" );
        }
        setters += ( i<multiplicity-1 ? ", " : ") {\n" );
      }
      for (int i=0; i<multiplicity; i++) {
        for (int j=0; j<_dimensions; j++) {
          setters += indent +"  "+"_F["+j+"]["+(offset+i)+"]=v"+i+j+";\n";
        }
      }
      setters += indent +"}\n\n";
    }
    return setters;
  }
  
  private String createFluxesSetters() {
    String setters = "";
    
    int offset = 0;
    for (String identifier : _variablesMap.keySet()) {
      int multiplicity = _variablesMap.get(identifier);
      setters = appendFluxesSetters(setters, identifier, multiplicity, offset);
      offset += multiplicity;
    }
    
    return setters;
  }

  public void writeHeader(BufferedWriter writer, String solverName, String projectName)
      throws java.io.IOException {
    String content = IOUtils.convertRessourceContentToString(
        "eu/exahype/variables/templates/VariablesHeader.template");
    
    content = content.replaceAll("\\{\\{Project\\}\\}", projectName);
    content = content.replaceAll("\\{\\{Solver\\}\\}", solverName);

    String variablesGetters = createVariablesGetters();
    String variablesSetters = createVariablesSetters();
    content = content.replaceAll("\\{\\{VariablesGetters\\}\\}", variablesGetters);
    content = content.replaceAll("\\{\\{VariablesSetters\\}\\}", variablesSetters);
    
    String fluxesGetters = createFluxesGetters();
    String fluxesSetters = createFluxesSetters();
    content = content.replaceAll("\\{\\{FluxesGetters\\}\\}", fluxesGetters);
    content = content.replaceAll("\\{\\{FluxesSetters\\}\\}", fluxesSetters);

    writer.write(content);
    writer.flush();
  }
}
