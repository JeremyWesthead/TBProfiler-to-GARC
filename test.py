import re

del_ = re.compile(r"""
                n|c\. #Leading type
                (-?[0-9]+) #Start pos
                (_-?[0-9]+)? #Maybe end pos
                del #Del
                ([AGCT]+)? #Bases deleted
                """, re.VERBOSE)
mut = "c.106_278del"
