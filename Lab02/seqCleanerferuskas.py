class DNAstring(str):
    # Define a new class that extends the built-in str class.
    # This new class is called DNAstring.
    
    def length(self):
        # Define a new method called length that returns the length of the string.
        # This method is a member of the DNAstring class, and takes the 'self' parameter.
        
        return len(self)
        # Return the length of the string using the built-in len function.
    
    def purify(self):
        # Define a new method called purify that returns a cleaned up version of the string.
        # This method is a member of the DNAstring class, and takes the 'self' parameter.
        
        purified = ''
        # Create an empty string to hold the purified version of the input string.
        
        count = 0
        # Create a counter to keep track of the number of 'N's that occur in a row.
        
        prev_char = ''
        # Create a variable to store the previous character in the input string.
        # This is used to detect runs of 'N's that need to be replaced with a count.
        
        for char in self:
            # Loop over each character in the input string.
            
            if char.islower():
                # If the character is lowercase, convert it to uppercase.
                char = char.upper()
                
            if char == 'N':
                # If the character is 'N', increment the count.
                count += 1
                
            else:
                if prev_char == 'N':
                    # If the previous character was 'N', add a count of 'N's to the purified string.
                    purified += '{%d}' % count
                    # Use string formatting to add the count to the purified string.
                    count = 0
                    # Reset the count to zero.
                    
                purified += char
                # Add the current character to the purified string.
                
            prev_char = char
            # Store the current character as the previous character for the next iteration.
            
        if prev_char == 'N':
            # If the last character in the string was 'N', add a count of 'N's to the purified string.
            purified += '{%d}' % count
            
        return purified
        # Return the purified string.
        
        
def main():
    ''' Get user DNA data and clean it up.'''
    # Define a function called main that prompts the user for DNA data, cleans it up, and prints the result.
    
    data = input('DNA data?')
    # Prompt the user for input.
    
    while data:
        # Loop until the user enters an empty string.
        
        thisDNA = DNAstring(data)
        # Create a new instance of the DNAstring class using the user's input string.
        
        pureData = thisDNA.purify()
        # Clean up the input string using the purify method of the DNAstring instance.
        
        print(pureData)
        # Print the purified string to the console.
        
        data = input('DNA data?')
        # Prompt the user for input again.
        
        
main()
# Call the main function to start the program. 