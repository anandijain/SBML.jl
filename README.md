# SBML

## http://co.mbine.org/standards/sbml

```xml
<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
   <model id="MyModel">
      <listOfFunctionDefinitions>
         zero or more
         <functionDefinition>...</functionDefinition>
         elements
      </listOfFunctionDefinitions>
      <listOfUnitDefinitions>
         zero or more
         <unitDefinition>...</unitDefinition>
         elements
      </listOfUnitDefinitions>
      <listOfCompartments>
         zero or more
         <compartment>...</compartment>
         elements
      </listOfCompartments>
      <listOfSpecies>
         zero or more
         <species>...</species>
         elements
      </listOfSpecies>
      <listOfParameters>
         zero or more
         <parameter>...</parameter>
         elements
      </listOfParameters>
      <listOfInitialAssignments>
         zero or more
         <initialAssignment>...</initialAssignment>
         elements
      </listOfInitialAssignments>
      <listOfRules>zero or more elements of subclasses ofRule</listOfRules>
      <listOfConstraints>
         zero or more
         <constraint>...</constraint>
         elements
      </listOfConstraints>
      <listOfReactions>
         zero or more
         <reaction>...</reaction>
         elements
      </listOfReactions>
      <listOfEvents>
         zero or more
         <event>...</event>
         elements
      </listOfEvents>
   </model>
</sbml>
```