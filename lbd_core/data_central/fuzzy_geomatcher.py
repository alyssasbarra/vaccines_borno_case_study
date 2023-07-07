import pandas as pd
import os, sys
import pyarabic
import re
from fuzzywuzzy import process

def placename_matcher(file1,
                      file2,
                      top_level_col1,
                      top_level_col2,
                      output_file = None,
                      second_level_col1 = None,
                      second_level_col2 = None,
                      third_level_col1 = None,
                      third_level_col2 = None,
                      header1 = 0,
                      header2 = 0,
                     ):

    '''
    file1 should be the vector or file column you are matching FROM, and file2 is the vector/column you are matching TO.
    In other words, file2 is the thing you are trying to get matches for. The program will look for the closest matches
    to the names in file2, in file1, and then pull the file1 names into a new column of file2 called "fuzzy_match_top" when a best match is found.

    The code refers to the 3 levels of fuzzy matching as "district", "subdistrict", and "village" because those are the levels it was
    originally written for, but any kind of multi-level data can be used. For example one could match top level=year, second level=sex and third level=admin2

    "District" is meant to be a category within which the fuzzy matcher will look for matches. For example,
    A user could specify Washington, and the program would only look for matches within Washington states for
    County matches, rather than comparing to every county in the US. This cuts down on computational time significantly.
    One could also specify a country to look for some lower admin unit in.

    Multilevel matching is implemented to allow you to only fuzzy match within a subset of data, which can significantly reduce computational time.
    For example, if you have a lot of admin2's, you could first fuzzy match on the admin1 level for top level fuzzy match, and then the program will use that
    closest match to only fuzzy match admin2's within an admin1. So instead of taking a long time to check an item in your list against EVERY item in the database
    you are checking it with, it will only check admin2's that have the same admin1. This can significantly reduce the likelihood that an incorrect fuzzy match will be found
    in a completely different area of a country.

    I'd recommend that you still check over the matches by sorting least close -> closest match for the whole dataset, to make sure that your non-100% matches are getting the right thing.
    Basically, this is a better vlookup.

    '''

        #reading in the data
    df1 = file1.copy()
    df2 = file2.copy()

    #create fuzzy matching function
    def fuzzy_match2(x, match_set):
        match, score = process.extractOne(x, match_set)
        return (match, score)

    df2[top_level_col2] = df2[top_level_col2].map(str)
    df1[top_level_col1] = df1[top_level_col1].map(str)


    districts = list(df2[top_level_col2].unique())

    #if user includes third_level
    if third_level_col1:
        df1[third_level_col1] = df1[third_level_col1].fillna('')

    fuzz_set1 = set(list(df1[top_level_col1]))

    for district in districts:
        district_subset = df2[df2[top_level_col2] == str(district)]
        subset_indices = district_subset.index.tolist()
        district_match = fuzzy_match2(district, fuzz_set1)
        df2.loc[subset_indices,'fuzzy_match_top'] = district_match[0]
        df2.loc[subset_indices,'score_top'] = district_match[1]

    if second_level_col1:
        #create a set of all the datapoints with a given district number
        for district in districts:
            fuzz_set = set(list(df1.loc[df1[top_level_col1] == district][second_level_col1]))

            #mark missing districts in "issue_level"
            if len(fuzz_set) == 0:
                df2.loc[df2[top_level_col2] == district, 'issue_level'] = 'DIST'
                print('No matches were found for district: {}'.format(district))
                continue

            #do the actual fuzzy match, on a given district's subdistrict.
            mask = df2[top_level_col2] == district
            df2.loc[mask, 'result'] = df2.loc[mask, second_level_col2].apply(lambda x: fuzzy_match2(x, fuzz_set))
            df2[['fuzzy_match_subdist', 'score_subdist']] = pd.DataFrame(df2['result'].tolist(),index=df2.index)

        subdistricts = df2[df2['score_subdist'] >= 0].fuzzy_match_subdist.unique().tolist()

        if third_level_col1:
            #now do matches on all city/village-level
            for subdistrict in subdistricts:

                # creates a set of the subdistricts which are in the fuzzy_match_subdist list.
                fuzz_set_sub = set(list(df1.loc[df1[second_level_col1] == subdistrict][third_level_col1]))

                mask2 = df2['fuzzy_match_subdist'] == subdistrict
                df2.loc[mask2, 'result2'] = df2.loc[mask2, third_level_col2].apply(lambda x: fuzzy_match2(x, fuzz_set_sub))
                df2[['fuzzy_match_village', 'score_village']] = pd.DataFrame(df2['result2'].tolist(),index=df2.index)
                #print("Villages in Subdistrict {} has been matched".format(subdistrict))

    #perfect_matches =  df2[['score_subdist'==100 & 'score_village'==100]]
    if output_file == None:
        return df2
    else:
        df2.to_excel(output_file, index=False)
