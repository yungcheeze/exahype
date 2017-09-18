package minitemp.syntaxtree;

import java.util.List;
import minitemp.Context;

/**
 * Abstract class to build a syntax tree from the token list.
 * The template can then be rendered by the syntax tree performing a DFS
 */
public abstract class SyntaxTree {
  
  /** parent node if available */
  public SyntaxTree parent;
  protected List<SyntaxTree> children = null;
  
  /** Add a node to this node */
  public void addNode(SyntaxTree node) {}
  
  /** render this node */
  public String render(Context context) throws IllegalArgumentException { return "";}
  
  /**
   * Render the children of this node.
   * It is expected that the node classes will use some kind of DFS when performing render
   *
   * @param the context of the rendering used to evaluate nodes and leafs
   * @return the rendered syntax tree
   */
  public String renderChildren(Context context) throws IllegalArgumentException {
    String result = "";
    for(SyntaxTree child : children) {
      result += child.render(context);
    }
    return result;
  }
  
  /**
   * Build a SyntaxTree from a template (String) and a regex to tokenize the template
   * Use the static Token.tokenize to tokenize the template with the regex
   *
   * @param the template to transform into a SyntaxTree
   * @param the regex to tokenize the template
   * @return a SyntaxTree ready to be rendered
   */
  static public final SyntaxTree buildTree(String template, String regex) 
      throws IllegalArgumentException
  {
    Token[] tokens = Token.tokenize(template, regex);
    
    SyntaxTree currentNode = new RootNode(); //start with a root node
    for(Token t : tokens) { //build the tree token by token
      if(t.type == Token.VAR_TOKEN) {             //variable is a Leaf
        currentNode.addNode(new VariableLeaf(t)); 
      } else if(t.type == Token.TEXT_TOKEN) {     //text is a Leaf
        currentNode.addNode(new TextLeaf(t));     
      } else if (t.type == Token.IF_OPEN_TOKEN) { //if logic start a subtree
        SyntaxTree newNode = new IfNode(t, currentNode);
        currentNode.addNode(newNode);
        currentNode = newNode;
      } else if (t.type == Token.IF_ELSE_TOKEN) { //else stay in the if subtree but notify the node
        if (currentNode instanceof IfNode) { // the current node has to be an IfNode
          IfNode p = (IfNode)(currentNode); 
          p.startElseBlock();
        } else {
          throw new IllegalArgumentException("Read an else block without being in an if block");
        }
      } else if (t.type == Token.IF_CLOSE_TOKEN) { //endif, go back to the parent of the subtree
        if(currentNode.parent == null) {
          throw new IllegalArgumentException("Read an endif block without being in an if block");
        }
        currentNode = currentNode.parent;
      }
    }
    
    return currentNode;
  }
  
}
