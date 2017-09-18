package minitemp.syntaxtree;

import minitemp.TemplateEngine;

public class Token {
  
  public static final int VAR_TOKEN      = 0; /** type of variable grammar token*/
  public static final int TEXT_TOKEN     = 1; /** type of text token*/
  public static final int IF_OPEN_TOKEN  = 2; /** type of logic grammar token: if */
  public static final int IF_ELSE_TOKEN  = 3; /** type of logic grammar token: else */
  public static final int IF_CLOSE_TOKEN = 4; /** type of logic grammar token: endif */
  
  
  public int type = -1; /** type of token, matched with the constant, by default -1 = invalid */
  private String rawContent; /** the token full string */
  
  
  public Token(String tokenAsString) {
    this.rawContent = tokenAsString;
    defineType();
  }
  
  
  private void defineType() {
    if(rawContent.startsWith(TemplateEngine.VAR_TOKEN_START)) {
      type = VAR_TOKEN;
    } else if(rawContent.startsWith(TemplateEngine.BLOCK_TOKEN_START)) {
      String tag = rawContent;
      if(tag.startsWith(TemplateEngine.STRIP_BLOCK_TOKEN_START)) {
        tag = tag.substring(TemplateEngine.STRIP_BLOCK_TOKEN_START.length()).trim();
      } else {
        tag = tag.substring(TemplateEngine.BLOCK_TOKEN_START.length()).trim();
      }
      if(tag.startsWith(TemplateEngine.LOGIC_ENDIF_TAG)) {
        type = IF_CLOSE_TOKEN;
      } else if(tag.startsWith(TemplateEngine.LOGIC_ELSE_TAG)) {
        type = IF_ELSE_TOKEN;
      } else if(tag.startsWith(TemplateEngine.LOGIC_IF_TAG)) {
        type = IF_OPEN_TOKEN;
      }
    } else {
      type = TEXT_TOKEN;
    }
  }
  
  @Override
  public String toString() {
    return type+" | \""+rawContent+"\"";
  }
  
  /**
   * Return the clean content of a token
   * Text token: everything in the token
   * Grammar token: the inside of the delimiters - tag + trim()
   */
  public String getContentClean() {
    if(type == VAR_TOKEN || type == IF_ELSE_TOKEN || type == IF_CLOSE_TOKEN) {
      String clean = rawContent.substring(TemplateEngine.BLOCK_TOKEN_START.length());
      clean = clean.substring(0, clean.length()-TemplateEngine.BLOCK_TOKEN_END.length());
      return clean.trim();
    } else if (type == TEXT_TOKEN) {
      return rawContent;
    } else if(type == IF_OPEN_TOKEN) {
      String clean = rawContent;
      if(clean.startsWith(TemplateEngine.STRIP_BLOCK_TOKEN_START)) 
      {
        clean = clean.substring(TemplateEngine.STRIP_BLOCK_TOKEN_START.length());
      } else {
        clean = clean.substring(TemplateEngine.BLOCK_TOKEN_START.length());
      }
      clean = clean.substring(0, clean.length()-TemplateEngine.BLOCK_TOKEN_END.length());
      clean = clean.trim().substring(TemplateEngine.LOGIC_IF_TAG.length()).trim(); //remove the tag
      return clean;
    }
    
    return "";
  }
  
  /**
   * Breaks a string into tokens using a given regex.
   * Apply strip token logic before building the token so return value != input splitted
   *
   * @param the template to evaluate as String
   * @param the regex used to split the template into tokens
   * @return the string splitted into token
   */
  static public Token[] tokenize(String template, String regex){       
    String[] tokens_s = template.split(regex);
    
    //apply strip block
    for(int i=0; i<tokens_s.length; i++) {
      if(tokens_s[i].startsWith(TemplateEngine.STRIP_BLOCK_TOKEN_START)) {
        if(i>0) {
          //remove trailing vertical whitespace from previous token
          tokens_s[i-1] = tokens_s[i-1].replaceAll("\\h*$", ""); 
        }
        if(i<tokens_s.length-1) {
          //remove leading whitespace + newline from next token
          tokens_s[i+1] = tokens_s[i+1].replaceAll("^\\h*\\R?", ""); 
        }
      }
    }
    
    //build token
    Token[] tokens = new Token[tokens_s.length];
    for(int i=0; i<tokens_s.length; i++) {
      tokens[i] = new Token(tokens_s[i]);
    }
    
    return tokens;
  }
  
}
