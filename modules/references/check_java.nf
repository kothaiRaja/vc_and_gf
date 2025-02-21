process CHECK_JAVA {
    tag "Check Java"
    container null  
	
    output:
    path "java_check.log", emit: java_output
	
    when:
    true  

    script:
    """
    # Get the Java version
    java_version=\$(java -version 2>&1 | grep -oP '(?<=version ")([0-9]+)' | head -1)
    
    if [[ -z "\$java_version" ]]; then
        echo " ERROR: Java is not installed or not in PATH."
        echo " ERROR: Java is not installed or not in PATH." > java_check.log
    elif [[ \$java_version -lt 21 ]]; then
        echo " WARNING: Java version \$java_version detected. Please update to version 21 or higher."
        echo " WARNING: Java version \$java_version detected. Please update to version 21 or higher." > java_check.log
    else
        echo " Java version \$java_version detected. Java is up-to-date."
        echo " Java version \$java_version detected. Java is up-to-date." > java_check.log
    fi
    """
}