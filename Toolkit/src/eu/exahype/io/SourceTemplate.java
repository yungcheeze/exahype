package eu.exahype.io;

import eu.exahype.io.IOUtils;

/**
 * The SourceTemplate class supersedes IOUtils and hides a template engine.
 *
 *
 **/
abstract public class SourceTemplate {
  abstract public void put(String key, String value);
  abstract public String toString();
  
  public static SourceTemplate getDefaultImplementation(String inline_template) {
    // here we choose about the default implementation
    return new TextReplaceTemplate(inline_template);
  }

  public static SourceTemplate fromRessourceContent(String relativePath) {
    String template = IOUtils.convertRessourceContentToString(relativePath);
    return getDefaultImplementation(template);
  }
}
