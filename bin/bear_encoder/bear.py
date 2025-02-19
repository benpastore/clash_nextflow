
import subprocess
import os

# o : residues in loop 
# d/f/c/a : base pairing 
# :/m : mismatch

"""
LOOP	j	k	l	m	n	o	p	q	r	s	t	u	v	w	x	y	z	^
STEM	a	b	c	d	e	f	g	h	i	=								
STEM_branch	A	B	C	D	E	F	G	H	I	J								
LEFTINTERNALLOOP	?	!	"	#	$	%	&	'	(	)	+							
BULGELEFT	[																	
LEFTINTERNALLOOP_branch	?	K	L	M	N	O	P	Q	R	S	T	U	V	W				
BULGELFETBRANCH	{																	
RIGHTINTERNALLOOP	?	2	3	4	5	6	7	8	9	0	>							
BULGERIGHT	]																	
RIGHTINTERNALLOOP_branch	?	Y	Z	~	?	_	|	/	\	@								
BULGERIGTHBRANCH	}																	

"""

script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)

if not os.path.exists(f"{script_dir}/tmp_files") : 
    os.mkdir(f"{script_dir}/tmp_files")

def annotate_bp_bear_encoder(sequence, basepairing) : 

    basepairing_bear = basepairing.replace("&", "")
    sequence_bear = sequence.replace("&", "")
    my_input = f">A\n{sequence_bear}\n{basepairing_bear}"
    tmp = open(f"{script_dir}/tmp_files/{sequence}.fa", 'w')
    tmp.write(my_input)
    tmp.close()


    result = subprocess.run(
        ["java", "-jar", f"{script_dir}/jar/BEAR_encoder.jar", "-i", f"{script_dir}/tmp_files/{sequence}.fa" ],
        text=True,          # Handle input/output as text.
        capture_output=True
    )
    output = result.stdout.strip().split("\n")
    bear_encoding = output[3]
    target_len = len(sequence.split("&")[0])

    bear_encoding_split = bear_encoding[0:target_len] + "&" + bear_encoding[target_len:len(sequence)]

    return bear_encoding_split
