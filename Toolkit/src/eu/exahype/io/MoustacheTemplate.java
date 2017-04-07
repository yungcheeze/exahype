package eu.exahype.io;

import eu.exahype.io.SourceTemplate;

import com.github.mustachejava.Mustache;
import com.github.mustachejava.DefaultMustacheFactory;
import com.github.mustachejava.MustacheFactory;

import java.io.StringReader;
import java.io.StringWriter;
import java.util.HashMap;
import java.util.Map;


/**
 * A template instance using Moustache-java,
 *
 *
 **/
public class MoustacheTemplate extends SourceTemplate {
  Mustache mustache;
  Map<String,Object> values;

  public MoustacheTemplate(String inline_template) {
     values = new HashMap<String,Object>();
     MustacheFactory mf = new DefaultMustacheFactory();
     mustache = mf.compile(new StringReader(inline_template), "example");
  }

  public void put(String key, String value) {
    values.put(key, value);
  }

  public String toString() {
    StringWriter sw = new StringWriter();
    mustache.execute(sw, values);
    return sw.toString();
  }
}
