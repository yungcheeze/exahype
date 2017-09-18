package minitemp.syntaxtree;

import minitemp.Context;

/**
 * Syntax tree leaf for variable '{{ var }}'
 * Currently only work with simple string replacement
 */
public class VariableLeaf extends SyntaxTree {
  
  /** Store the token for later rendering */
  private Token token;
  
  public VariableLeaf(Token t) {
    token = t;
  }

  @Override
  public String render(Context c) throws IllegalArgumentException {    
    return c.evaluateString(token.getContentClean());
  }
  
}
