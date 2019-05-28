<?php

	//echo "Your job is being processed with specified job id. <br/> ";
	echo "";
	

if (isset($argc)) {

	echo "The evaluation result is available <a href = \"". $argv[1]. "\">here. </a> <br/>" ;

	echo "Log file is <a href =\"". $argv[2] . "\">here. </a> and the job directory is <a href =\"". $argv[3] . "\">here. </a> <br/> " ;
	
	echo "Usually it takes  a few seconds(depends on the sequence length).<br/>"	;
}
else {
	echo "argc and argv disabled\n";
}

?>