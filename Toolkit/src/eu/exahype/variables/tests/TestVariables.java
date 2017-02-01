package eu.exahype.variables.tests;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.LinkedHashMap;
import java.util.Map;

import eu.exahype.variables.Variables;

public class TestVariables {
  
  public static void testVariables() {
    int dimensions=2;
    
    Map<String,Integer> variablesMap = new LinkedHashMap<String,Integer>();
    variablesMap.put("rho", 1);
    variablesMap.put("u", 3);
    variablesMap.put("E", 1);
    
    Map<String,Integer> parametersMap = new LinkedHashMap<String,Integer>();
    parametersMap.put("matScalar", 1);
    parametersMap.put("matVector", 3);
    
    Variables variables = new Variables(variablesMap,parametersMap,dimensions);
    
    BufferedWriter bufferedWriter = new BufferedWriter(new OutputStreamWriter(System.out));
    
    try {
      System.out.println("Content of header file: ");
      variables.writeHeader(bufferedWriter, "MySolver", "MyProject");
      bufferedWriter.flush();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
  
  public static void main(String[] args) {
    testVariables();
  }
}
