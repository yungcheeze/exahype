package minitemp.syntaxtree;

import minitemp.Context;

/**
 * Syntax tree node for Branch logic.
 * Store what's contained between its branches in subtree and render the correct one
 */
public class IfNode extends SyntaxTree {
  
  private Token condition;       //condition for the if
  private RootNode ifBlock;      //subtree for the if block
  private RootNode elseBlock;    //subtree for the else block
  private RootNode currentBlock; //pointer to the subtree being filled
  
  public IfNode(Token t, SyntaxTree p) {
    condition = t;
    this.parent = p;
    ifBlock = new RootNode();
    elseBlock = new RootNode();
    currentBlock = ifBlock;
  }
  
  /** Add a node to the current subtree */
  @Override
  public void addNode(SyntaxTree t) {
    this.currentBlock.addNode(t);
  }
  
  /** Trigger to start filling the else subtree, triggered when encountering the relevant else */
  public void startElseBlock() {
    currentBlock = elseBlock;
  }
  
  //TODO extend to boolean expression
  /** Render the subtree corresponding to the correct branch (eventually empty) */
  @Override
  public String render(Context c) throws IllegalArgumentException {
    boolean value = c.evaluateBoolean(condition.getContentClean());
    if(value) {
      return ifBlock.render(c);
    } else {
      return elseBlock.render(c);
    }
  }
  
}
