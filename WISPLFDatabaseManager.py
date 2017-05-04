import datetime
import numpy as np
import sqlite3


class WISPLFDatabaseManager:
    """
    Utility class to manage persistence of line identification data during the
    WISP line-finding process.
    """

    validFlags = [
        'REJECT',
        'CONTAM',
        'ZEROTH',
        'CONTIN',
        'MISC'
    ]

    fitResultKeys = ['redshift',
                     'redshift_err',
                     'dz_oiii',
                     'dz_oii',
                     'dz_siii_he1',
                     'fwhm_g141',
                     'fwhm_g141_err',
                     'oii_flux',
                     'oii_error',
                     'oii_ew_obs',
                     'hg_flux',
                     'hg_error',
                     'hg_ew_obs',
                     'hb_flux',
                     'hb_error',
                     'hb_ew_obs',
                     'oiii_flux',
                     'oiii_error',
                     'oiii_ew_obs',
                     'hanii_flux',
                     'hanii_error',
                     'hanii_ew_obs',
                     'sii_flux',
                     'sii_error',
                     'sii_ew_obs',
                     'siii_9069_flux',
                     'siii_9069_error',
                     'siii_9069_ew_obs',
                     'siii_9532_flux',
                     'siii_9532_error',
                     'siii_9532_ew_obs',
                     'he1_flux',
                     'he1_error',
                     'he1_ew_obs']

    def __init__(self, dbFileNamePrefix):
#        print('Using sqlite3 version {}'.format(sqlite3.version))
        self.dbFileNamePrefix = dbFileNamePrefix
        self.dbFileName = '{}_sqlite.db'.format(self.dbFileNamePrefix)
        self.dbConnection = sqlite3.connect(self.dbFileName)
        self.dbConnection.row_factory = sqlite3.Row
        self.dbCursor = self.dbConnection.cursor()
        self.checkAndInitTables()

    def __del__(self):
        self.dbConnection.commit()
        self.dbConnection.close()

    def checkAndInitTables(self):
        self.createCatalogueTable()
        self.createAnnotationTable()
        self.createFlagTable()

    def createCatalogueTable(self):
#        print('Creating catalogue table in SQLLite database...')
        self.dbCursor.execute('''CREATE TABLE IF NOT EXISTS catalogue (
                              ParNo int,
                              ObjID int,
                              RA real,
                              Dec real,
                              Jmagnitude real,
                              Hmagnitude real,
                              A_IMAGE real,
                              B_IMAGE real,
                              redshift real,
                              redshift_err real,
                              dz_oiii real,
                              dz_oii real,
                              dz_siii_he1 real,
                              fwhm_g141 real,
                              fwhm_g141_err real,
                              oii_flux real,
                              oii_error real,
                              oii_ew_obs real,
                              hg_flux real,
                              hg_error real,
                              hg_ew_obs real,
                              hb_flux real,
                              hb_error real,
                              hb_ew_obs real,
                              oiii_flux real,
                              oiii_error real,
                              oiii_ew_obs real,
                              hanii_flux real,
                              hanii_error real,
                              hanii_ew_obs real,
                              sii_flux real,
                              sii_error real,
                              sii_ew_obs real,
                              siii_9069_flux real,
                              siii_9069_error real,
                              siii_9069_ew_obs real,
                              siii_9532_flux real,
                              siii_9532_error real,
                              siii_9532_ew_obs real,
                              he1_flux real,
                              he1_error real,
                              he1_ew_obs real,
                              ContamFlag int,
                              EntryTime text
                              )''')
        self.dbConnection.commit()
#        print('Done.' if self.dbCursor.rowcount >
#              0 else 'Table was already present.')

    def createAnnotationTable(self):
#        print('Creating annotations table in SQLLite database...')
        self.dbCursor.execute('''CREATE TABLE IF NOT EXISTS annotations (
                         ParNo int,
                         ObjID int,
                         Comment text
                         )''')
        self.dbConnection.commit()
#        print('Done.' if self.dbCursor.rowcount >
#              0 else 'Table was already present.')

    def createFlagTable(self):
#        print('Creating annotations table in SQLLite database...')
        self.dbCursor.execute('''CREATE TABLE IF NOT EXISTS flags (
                         ParNo int,
                         ObjID int,
                         FlagName text,
                         FlagValue int
                         )''')
        self.dbConnection.commit()
#        print('Done.' if self.dbCursor.rowcount >
#              0 else 'Table was already present.')

    def saveCatalogueEntry(self, catalogueEntryData):
        query = 'INSERT INTO catalogue VALUES ({})'.format(
            ','.join(['?'] * len(catalogueEntryData))
        )
        # print(query)
        self.dbCursor.execute(query, catalogueEntryData)
        self.dbConnection.commit()

    def loadCatalogueEntry(self, parNumber, objectId):
        query = 'SELECT * FROM catalogue WHERE (ParNo=? AND ObjID=?)'
        self.dbCursor.execute(query, (parNumber, objectId))
        catalogueEntryData = self.dbCursor.fetchall()
        if len(catalogueEntryData) < 1:
            return None
        # FIXME: For now, just return the last row added
        nonFitResults = [catalogueEntryData[-1][key]
                         for key in catalogueEntryData[-1].keys() if key not in WISPLFDatabaseManager.fitResultKeys]
        fitResults = {key: catalogueEntryData[-1][key]
                      for key in WISPLFDatabaseManager.fitResultKeys}
        return tuple(nonFitResults), fitResults

    def getMostRecentObject(self, parNumber=None):
        query = '''SELECT ParNo, ObjID, EntryTime
        FROM catalogue{}
        ORDER BY DATETIME(EntryTime)
        DESC LIMIT 1'''.format('WHERE ParNo = ?' if parNumber is not None else '')
        self.dbCursor.execute(
            *[argument for argument in [query, parNumber] if argument is not None])
        mostRecentEntry = self.dbCursor.fetchone()
        if mostRecentEntry is not None:
            return mostRecentEntry['ParNo'], mostRecentEntry['ObjId'], mostRecentEntry['EntryTime']
        return None

    def layoutCatalogueData(self,
                            parNumber,
                            objectId,
                            ra,
                            dec,
                            jMagnitude,
                            hMagnitude,
                            aImage,
                            bImage,
                            fitResults,
                            flagContent):
        if fitResults is None:
            fitResults = {
                key: None for key in WISPLFDatabaseManager.fitResultKeys}
        return (parNumber, objectId, ra, dec, jMagnitude, hMagnitude, aImage, bImage,
                fitResults['redshift'],
                fitResults['redshift_err'],
                fitResults['dz_oiii'],
                fitResults['dz_oii'],
                fitResults['dz_siii_he1'],
                fitResults['fwhm_g141'],
                fitResults['fwhm_g141_err'],
                fitResults['oii_flux'],
                fitResults['oii_error'],
                fitResults['oii_ew_obs'],
                fitResults['hg_flux'],
                fitResults['hg_error'],
                fitResults['hg_ew_obs'],
                fitResults['hb_flux'],
                fitResults['hb_error'],
                fitResults['hb_ew_obs'],
                fitResults['oiii_flux'],
                fitResults['oiii_error'],
                fitResults['oiii_ew_obs'],
                fitResults['hanii_flux'],
                fitResults['hanii_error'],
                fitResults['hanii_ew_obs'],
                fitResults['sii_flux'],
                fitResults['sii_error'],
                fitResults['sii_ew_obs'],
                fitResults['siii_9069_flux'],
                fitResults['siii_9069_error'],
                fitResults['siii_9069_ew_obs'],
                fitResults['siii_9532_flux'],
                fitResults['siii_9532_error'],
                fitResults['siii_9532_ew_obs'],
                fitResults['he1_flux'],
                fitResults['he1_error'],
                fitResults['he1_ew_obs'],
                flagContent,
                str(datetime.datetime.now().isoformat())
                )

    def saveAnnotation(self, annotationData):
        query = 'INSERT INTO annotations VALUES ({})'.format(
            ','.join(['?'] * len(annotationData))
        )
        self.dbCursor.execute(query, annotationData)
        self.dbConnection.commit()

    def setFlagsFromString(self, parNumber, objId, flagDataString, delimiter=','):
        # Assumes that string is a delimiter separated list of flags and values
        flagDataTokens = [token.strip()
                          for token in flagDataString.split(delimiter)]
        flagNames = flagDataTokens[::2]
        flagValues = flagDataTokens[1::2]
        flagData = list(zip(flagNames, flagValues))
        self.setFlags(parNumber, objId, flagData)

    def setFlags(self, parNumber, objId, flagData):
        # Only attempt to set values for valid flags
        flagData = [(parNumber, objId, flagDatum[0], flagDatum[1])
                    for flagDatum in flagData
                    if flagDatum[0] in WISPLFDatabaseManager.validFlags]
        # flagData can be a list of (flagName, flagValue) tuples
        query = 'INSERT INTO flags VALUES (?, ?, ?, ?)'
        self.dbCursor.executemany(query, flagData)
        self.dbConnection.commit()

    def writeCatalogueTextFile(self):
        pass
