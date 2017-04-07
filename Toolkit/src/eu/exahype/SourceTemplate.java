package eu.exahype;

import eu.exahype.IOUtils;

/**
 * The SourceTemplate class supersedes IOUtils and hides a template engine.
 *
 * This class is introduced for the migration to a fully featured templating engine.
 * At first, it will just supersede what IOUtils did before.
 *
 *
 **/
public class SourceTemplate {
  String tpl;

  public SourceTemplate(String template) {
    // for the beginning, just masquerade what a String would do
    this.tpl = template;
  }

  public void put(String key, String value) {
    this.tpl = this.tpl.replaceAll("\\{\\{"+key+"\\}\\}", value);
  }

  public String toString() {
    return tpl;
  }

  public static SourceTemplate fromRessourceContent(String relativePath) {
    String template = IOUtils.convertRessourceContentToString(relativePath);
    return new SourceTemplate(template);
  }
}
