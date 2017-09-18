package minitemp.syntaxtree;

import java.util.ArrayList;
import minitemp.Context;

/**
 * Syntax tree node with multiple children. Do not contain a token.
 * Also the root of the syntax tree or any subtree
 */
public class RootNode extends SyntaxTree {
  
  public RootNode() {
    this.children = new ArrayList<SyntaxTree>();
  }
  
  @Override
  public void addNode(SyntaxTree t) {
    this.children.add(t);
  }
  
  /** Performs a DFS */
  @Override
  public String render(Context c) throws IllegalArgumentException {
    return renderChildren(c);
  }
  
}
