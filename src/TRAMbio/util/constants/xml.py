from enum import Enum
from typing import Dict, Tuple
import TRAMbio


################################
# Component XML Conversion #####
################################

_INDENT_WIDTH = 4

class XMLConstants(Enum):
    NODE_TAG = 'node'
    NODES_TAG = 'nodes'
    HALO_TAG = 'halo'
    COMPONENT_TAG = 'component'
    COMPONENTS_TAG = 'components'
    GRAPH_TAG = 'graph'
    STATE_TAG = 'state'
    STATES_TAG = 'states'

    STRUCTURE_ATTRIBUTE_NAME = 'structure'
    SIZE_ATTRIBUTE_NAME = 'size'
    KEY_ATTRIBUTE_NAME = 'key'

    K_ATTRIBUTE_NAME = 'k'
    L_ATTRIBUTE_NAME = 'l'
    CATEGORY_ATTRIBUTE_NAME = 'category'

    XML_INDENT = ' ' * _INDENT_WIDTH
    XML_INDENT2 = ' ' * (_INDENT_WIDTH * 2)
    XML_INDENT3 = ' ' * (_INDENT_WIDTH * 3)
    XML_INDENT4 = ' ' * (_INDENT_WIDTH * 4)
    XML_INDENT5 = ' ' * (_INDENT_WIDTH * 5)

    SPECIAL_CHARACTER = '#'

    XSD_COMPONENTS = f'xmlns="tram:components" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="tram:components https://raw.githubusercontent.com/gate-tec/TRAMbio/refs/heads/release/v{TRAMbio.__version__}/resources/components.xsd"'
    XSD_PEBBLE = f'xmlns="tram:pebble" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="tram:pebble https://raw.githubusercontent.com/gate-tec/TRAMbio/refs/heads/release/v{TRAMbio.__version__}/resources/pebble_game.xsd"'


XML_NAMESPACE = {'tram': 'tram:components'}


def components_namespace(const: XMLConstants) -> Tuple[str, Dict[str, str]]:
    return 'tram:' + const.value, XML_NAMESPACE


XML_SCHEME_COMPONENTS = r"""<?xml version="1.0" encoding="UTF-8"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema"
           targetNamespace="tram:components"
           xmlns="tram:components"
           xmlns:tram="tram:components"
           elementFormDefault="qualified"
           attributeFormDefault="unqualified">
  <!-- Type Definitions -->
  <xs:simpleType name="id.type">
    <xs:restriction base="xs:int">
      <xs:minInclusive value="1" />
    </xs:restriction>
  </xs:simpleType>
  <xs:simpleType name="key.type">
    <xs:restriction base="xs:float">
      <xs:pattern value="-INF|-?[0-9]+\.[0-9]+|0|[1-9][0-9]*" />
    </xs:restriction>
  </xs:simpleType>
  <xs:simpleType name="size.type">
    <xs:restriction base="xs:int">
      <xs:minInclusive value="0" />
    </xs:restriction>
  </xs:simpleType>
  <xs:simpleType name="node.type">
    <xs:restriction base="xs:string">
      <xs:pattern value="[A-Z][0-9]{4}.[ A-Z0-9]{2}[A-Z0-9]:[A-Z]+[A-Z0-9']{0,3}" />
    </xs:restriction>
  </xs:simpleType>
  <xs:simpleType name="structure.type">
    <xs:restriction base="xs:string">
      <xs:minLength value="1" />
    </xs:restriction>
  </xs:simpleType>
  <!-- Structure Definition -->
  <xs:element name="graph">
    <xs:complexType>
      <xs:sequence>
        <xs:element name="components">
          <xs:complexType>
            <xs:sequence>
              <xs:element name="component" maxOccurs="unbounded" minOccurs="1">
                <xs:complexType>
                  <xs:sequence>
                    <xs:element name="nodes">
                      <xs:complexType>
                        <xs:sequence>
                          <xs:element type="node.type" name="node" maxOccurs="unbounded" minOccurs="1"/>
                        </xs:sequence>
                      </xs:complexType>
                    </xs:element>
                    <xs:element name="halo">
                      <xs:complexType mixed="true">
                        <xs:sequence>
                          <xs:element type="node.type" name="node" maxOccurs="unbounded" minOccurs="0"/>
                        </xs:sequence>
                      </xs:complexType>
                    </xs:element>
                  </xs:sequence>
                  <xs:attribute type="id.type" name="id" use="required"/>
                  <xs:attribute type="size.type" name="size" use="required"/>
                  <xs:attribute type="structure.type" name="structure" use="optional"/>
                </xs:complexType>
              </xs:element>
            </xs:sequence>
            <xs:attribute type="size.type" name="size" use="required"/>
          </xs:complexType>
          <xs:unique name="components_component_id_unique">
            <xs:annotation>
              <xs:documentation xml:lang="en">
                Ensures: uniqueness of id attributes of &lt;component&gt; children of this &lt;components&gt; element.
              </xs:documentation>
            </xs:annotation>
            <xs:selector xpath="./tram:component"/>
            <xs:field xpath="@id"/>
          </xs:unique>
        </xs:element>
        <xs:element name="states">
          <xs:complexType>
            <xs:sequence>
              <xs:element name="state" maxOccurs="unbounded" minOccurs="1">
                <xs:complexType>
                  <xs:sequence>
                    <xs:element name="component" maxOccurs="unbounded" minOccurs="1">
                      <xs:complexType>
                        <xs:sequence>
                          <xs:element name="halo">
                            <xs:complexType>
                              <xs:sequence>
                                <xs:element type="node.type" name="node" maxOccurs="unbounded" minOccurs="0"/>
                              </xs:sequence>
                            </xs:complexType>
                          </xs:element>
                          <xs:element name="components">
                            <xs:complexType>
                              <xs:sequence>
                                <xs:element name="component" maxOccurs="unbounded" minOccurs="1">
                                  <xs:complexType>
                                    <xs:attribute type="id.type" name="id" use="required"/>
                                  </xs:complexType>
                                </xs:element>
                              </xs:sequence>
                            </xs:complexType>
                          </xs:element>
                        </xs:sequence>
                      </xs:complexType>
                    </xs:element>
                  </xs:sequence>
                  <xs:attribute type="key.type" name="key" use="required"/>
                </xs:complexType>
                <xs:unique name="state_component_id_unique">
                  <xs:annotation>
                    <xs:documentation xml:lang="en">
                      Ensures: uniqueness of id attributes of &lt;component&gt; children in a &lt;components&gt; list of any &lt;component&gt; child of this &lt;state&gt; element.
                    </xs:documentation>
                  </xs:annotation>
                  <xs:selector xpath=".//tram:components/tram:component"/>
                  <xs:field xpath="@id"/>
                </xs:unique>
              </xs:element>
            </xs:sequence>
            <xs:attribute type="key.type" name="key" use="optional"/>
          </xs:complexType>
          <xs:unique name="states_state_key_unique">
            <xs:annotation>
              <xs:documentation xml:lang="en">
                Ensures: uniqueness of key attributes of &lt;state&gt; of this &lt;states&gt; element.
              </xs:documentation>
            </xs:annotation>
            <xs:selector xpath="./tram:state"/>
            <xs:field xpath="@key"/>
          </xs:unique>
        </xs:element>
      </xs:sequence>
    </xs:complexType>
  </xs:element>
</xs:schema>"""
