default: fsps

fsps: force_look
	(cd fsps; make _fsps)

force_look:
	true

