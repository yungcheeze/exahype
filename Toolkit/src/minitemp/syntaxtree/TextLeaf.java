package minitemp.syntaxtree;

import minitemp.Context;

/**
 * Syntax tree leaf for pure text without any special meaning
 */
public class TextLeaf extends SyntaxTree {
  
  private String text;
  
  public TextLeaf(Token t) {
    text = t.getContentClean();
  }

  @Override
  public String render(Context c) throws IllegalArgumentException {
    return text;
  }
  
}
