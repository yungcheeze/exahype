package minitemp;

import minitemp.syntaxtree.Token;
import minitemp.syntaxtree.SyntaxTree;

/**
 * Main class
 *
 * Once the context is build and initialized, use render(template, context) to run the template
 * engine on the template with the given context
 *
 * You can change the grammar by editing the constant from this class
 */
public class TemplateEngine {
  
  //Configuration parameters
  //------------------------
  /** Var tokens contain variable and come alone */
  public static final String VAR_TOKEN_START   = "{{";
  public static final String VAR_TOKEN_END     = "}}";
  /** Block tokens contain logic and come in group, requiring a syntax tree to evaluate */
  public static final String BLOCK_TOKEN_START = "{%";
  public static final String BLOCK_TOKEN_END   = "%}";
  /** Special block token start delimiter that strip whitespaces around it */
  public static final String STRIP_BLOCK_TOKEN_START = "{%-"; //must start like BLOCK_TOKEN_START
  
  /** Tags of the logic block for Branches */
  public static final String LOGIC_IF_TAG    = "if";
  public static final String LOGIC_ELSE_TAG  = "else";
  public static final String LOGIC_ENDIF_TAG = "endif";
  
  /** Tags of the logic block for Loops */
  public static final String LOGIC_FOR_TAG     = "for";
  public static final String LOGIC_FOR_SET_TAG = "in"; // {% for value in collection %}
  public static final String LOGIC_ENDFOR_TAG  = "endfor";
  
  /** The size of the content of a grammar token has to be constrained for the look trick */
  private static final int TOKEN_MAX_SIZE = 99999;
  
  /** Regex used to tokenize the template */
  private String regex;
  
  public TemplateEngine() {
    buildRegex(); //set this.regex
  }
  
  public String getRegex() {
    return regex;
  }
  
  /**
   * Build a regex pattern to split a string (the template) into tokens (grammar and text ones)
   *
   * A grammar token is delimiter start, content, delimiter end
   * All the text between two grammar token (or before the first/ after the last) is one text token
   * Grammar tokens are contrained in size by TOKEN_MAX_SIZE (arbitrary large int)
   * Text tokens can be arbitrary long
   *
   * The base regex find the grammar tokens, the regex then use it with the lookahead and lookbehind 
   * trick to keep the grammar token matches in the split result
   */
  private void buildRegex() {
    //constrains the size of a grammar token for the look trick
    final String tokenContent = ".{0,"+TOKEN_MAX_SIZE+"}?";
    
    // tokens
    final String varToken = escRegex(VAR_TOKEN_START)+tokenContent+escRegex(VAR_TOKEN_END);
    final String blockToken = escRegex(BLOCK_TOKEN_START)+tokenContent+escRegex(BLOCK_TOKEN_END);
    
    //match one of the grammar tokens
    final String tokenPattern = varToken+"|"+blockToken; 
    
    //lookahead+lookbehind trick
    this.regex = "(?="+tokenPattern+")|(?<="+tokenPattern+")"; 
  }
  
  /** escape regex special char '{' and '}' by adding '\' */
  private String escRegex(String regex) {
    // ReplaceAll uses regex itself hence the '\\' and '\\\\'
    return regex.replaceAll("\\{", "\\\\{").replaceAll("\\}", "\\\\}");
  }  
  
  /**
   * Evaluate a template with a given context.
   *
   * @param the template to render as a String
   * @param the context to use to render the template (of Context class)
   * @return a string with the rendered template
   */
  public String render(String template, Context context) throws IllegalArgumentException {
    return SyntaxTree.buildTree(template, this.regex).render(context);
  }
  
}
