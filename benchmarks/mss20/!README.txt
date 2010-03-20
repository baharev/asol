To solve this problem with ASOL, you have to enable the problem specific
splitting rule and disable index sainity check, then compile the solver from the
source.

Comment out the following line in the main.cpp
    index_sanity_check(vars, index, n_vars, w);
	
and uncomment the following lines in the main.cpp

    if (!toDelete && index != -1 ) {
        select_custom_index(box_orig, n, w, index);
    }

Now, compile the solver.