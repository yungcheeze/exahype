# Minimalist Template Engine

A minimalist template engine in java with no dependency outside the standard java libraries (java >= 1.8).

/!\ Some inputs are not sanitize, use this template engine responsibly.

Github: https://github.com/gallardjm/minitemp

## Quick start

Hello World!

    String template = "Hello {{place}}{% if excited %}!{% if mad %}!!!1!!{% endif %}{% else %}.{% endif %}";

    TemplateEngine engine = new TemplateEngine();
    Context context = new Context();
    context.put("place", "world");
    context.put("excited", true);
    context.put("mad", false);

    String result = engine.render(template, context);

## Supported template syntax

The context contains the key->values pairs used by the expression in the grammar tokens, defined below. A context's key has to be a valid variable name (only characters from [a-zA-Z_0-9], at least one). Values can be of any java standard object type (e.g. String, Boolean, Double, Integer) or array of a standard object type. Note that lists are allowed as value and that primitive types will be casted to the corresponding object type automatically (e.g. int to Integer).

The Context uses a Javascript engine (from the package javax.script) to evaluate expression in the grammar tokens using the allocated key->value pairs. All expression inside the grammar tokens have to be valid javascript expressions that can use the keys as a variable with value the one associated to the key. 

/!\ Expressions and values are not sanitized.

### Text replacement token

Default delimiters: ```{{``` and ```}}```.

Any expression resulting in a String, an int or a double is valid inside these delimiters. The token is replaced by the output of the expression.

### Logic token

Default delimiters: ```{%``` and ```%}```.

These delimiters allow some logic in the template processing. Nested logic is supported.

The following logic is implemented:

* Branch ```{% if foo %} {% else %} {% endif %}```
  - foo is any boolean expression (javascript syntax)
  - else is optional
  
* Loop ```{% for foo in bar %}  {% endfor %}```
  - bar is any Collection of object
  - the content between the delimiters will be rendered for each different value of foo
  - if foo was already defined, it's value is overridden within the loop and restored afterward
  

### Misc

#### Strip logic block

With ```{%-``` as start delimiter of a logic token, strip all whitespaces before the token + all whitespaces and one newline if present after the token. This allow more readable templates without unecessary whitespaces and newlines in the result.

Example: ``` foo \n {%- x %}   \n bar``` is evaluated as ``` foo \n{% x %} bar```. 

