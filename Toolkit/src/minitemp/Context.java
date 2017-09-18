package minitemp;

import javax.script.ScriptEngineManager;
import javax.script.ScriptEngine;
import javax.script.ScriptException;

/**
 * Context class
 *
 * Store the pair of key->value and use it to evaluate expressions.
 * Type allowed for the value: String, boolean, int, double
 * A key can have only one value, redefining it (even with a value of a different type) override the
 * previous value
 *
 * Use put(key, value) to add a key and its value to the context.
 * evaluateX evaulate an expression and return a value of type X.
 */
public class Context {
  
  /** Use a Javascript engine to store the context and evaluate expressions */
  private ScriptEngine engine;
  /** 
   * If corrupted this context cannot perform any evaluation
   * Used to avoid throwing exceptions to handle during put operations
   */
  private boolean corrupted = false; //
  /** epsilon = 10^-11, use find int returned as double */
  private static final double EPS = 0.00000000001;
  
  public Context() {    
    engine = (new ScriptEngineManager()).getEngineByName("JavaScript");
  }
  
  /** 
   * Put a key->value pair in the context 
   * Errors are catched but trigger warnings and corrupt the Context (=unusable to render)
   */
  public void put(String key, Object value) {
    try {
      validateKey(key);
      engine.put(key, value); 
    } catch (Exception e) {
      System.err.println("ERROR in Context: put("+key+", "+value+") failed. Context is corrupted.");
      corrupted = true;
    }
  }
  
  /** Verify that a key is one word with normal characters ([a-zA-Z_0-9]+) */
  private void validateKey(String key) throws IllegalArgumentException {
    if(key == null || !key.matches("^\\w+$")) {
      System.err.println("ERROR in Context: invalid key: "+key);
      throw new IllegalArgumentException("Invalid key: "+key);
    }
  }
  
  /** Get a String from the context */
  public String evaluateString(String expression) throws IllegalArgumentException { 
    if(corrupted) {
      throw new IllegalArgumentException("Cannot evaluate corrupted context");
    }
    
    Object valueRaw;
    try {
      valueRaw = engine.eval(expression);
    } catch (ScriptException e) {
      throw new IllegalArgumentException("Cannot evaluate token's expression: "+expression);
    }
    
    if(valueRaw instanceof String) {
      return (String)valueRaw;
    } else if(valueRaw instanceof Double) {
      double value = (Double)valueRaw;
      String output = "";
      if(Math.abs(Math.floor(value)-value) < EPS) {
        output += ((int)value);
      } else {
        output = ""+value;
      }
      return output;
    }
    
    throw new IllegalArgumentException("Cannot evaluate token's expression as String: "+expression);
  }
  
  /** Get a boolean from the context */
  public boolean evaluateBoolean(String expression) throws IllegalArgumentException { 
    if(corrupted) {
      throw new IllegalArgumentException("Cannot evaluate corrupted context");
    }
    
    Object valueRaw;
    try {
      valueRaw = engine.eval(expression);
    } catch (ScriptException e) {
      throw new IllegalArgumentException("Cannot evaluate token's expression: "+expression);
    }
    
    if(valueRaw instanceof Boolean) {
      return ((Boolean)valueRaw).booleanValue();
    }
    
    throw new IllegalArgumentException("Cannot evaluate token's expression as boolean: "+expression);
  }
  
}
