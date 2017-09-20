package minitemp.syntaxtree;

import minitemp.Context;

/**
 * Syntax tree leaf for pure text without any special meaning
 */
public class TextLeaf extends SyntaxTree {
  
  private String text;
  
  public TextLeaf(Token token) {
    text = token.getContentClean();
  }

  @Override
  public String render(Context context) throws IllegalArgumentException {
    return text;
  }
  
}
