#!/usr/bin/env python3
"""rhea_text.reaction_title: substrate -> product from Rhea equation text."""
from rhea_text import reaction_title

assert reaction_title("N-acetylserotonin + S-adenosyl-L-methionine = melatonin + S-adenosyl-L-homocysteine + H(+)") == "N-acetylserotonin → melatonin"
assert reaction_title("melatonin + S-adenosyl-L-homocysteine + H(+) = N-acetylserotonin + S-adenosyl-L-methionine") == "N-acetylserotonin → melatonin"      # SAM side first
assert reaction_title("2 S-adenosyl-L-methionine + a thing = 2 S-adenosyl-L-homocysteine + a methylated thing + 2 H(+)") == "a thing → a methylated thing"
assert reaction_title("x + S-adenosyl-L-methionine => S-adenosyl-L-homocysteine + H(+)") is None                                    # nothing left on one side
assert len(reaction_title("a" * 80 + " + S-adenosyl-L-methionine = " + "b" * 80 + " + S-adenosyl-L-homocysteine")) <= 60
print("rhea_text OK")
