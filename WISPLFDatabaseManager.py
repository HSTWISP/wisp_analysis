import datetime
import numpy as np
import sqlite3


class WISPLFDatabaseManager:
    """
    Utility class to manage persistence of line identification data during the
    WISP line-finding process.
    """

    fitResultKeys = ['redshift',
                     'redshift_err',
                     'dz_oiii',
                     'dz_oii',
                     'dz_siii_he1',
                     'G141_FWHM_Obs',
                     'G141_FWHM_Obs_err',
                     'oii_flux',
                     'oii_error',
                     'oii_EW_obs',
                     'hg_flux',
                     'hg_err',
                     'hg_EW_obs',
                     'hb_flux',
                     'hb_err',
                     'hb_EW_obs',
                     'oiii_flux',
                     'oiii_err',
                     'oiii_EW_obs',
                     'hanii_flux',
                     'hanii_err',
                     'hanii_EW_obs',
                     'sii_flux',
                     'sii_err',
                     'sii_EW_obs',
                     'siii_9069_flux',
                     'siii_9069_err',
                     'siii_9069_EW_obs',
                     'siii_9532_flux',
                     'siii_9532_err',
                     'siii_9532_EW_obs',
                     'he1_10830_flux',
                     'he1_10830_err',
                     'he1_10830_EW_obs']

    def __init__(self, dbFileNamePrefix):
        print('Using sqlite3 version {}'.format(sqlite3.version))
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
        self.createCommentTable()

    def createCatalogueTable(self):
        print('Creating catalogue table in SQLLite database...')
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
                         G141_FWHM_Obs real,
                         G141_FWHM_Obs_err real,
                         oii_flux real,
                         oii_error real,
                         oii_EW_obs real,
                         hg_flux real,
                         hg_err real,
                         hg_EW_obs real,
                         hb_flux real,
                         hb_err real,
                         hb_EW_obs real,
                         oiii_flux real,
                         oiii_err real,
                         oiii_EW_obs real,
                         hanii_flux real,
                         hanii_err real,
                         hanii_EW_obs real,
                         sii_flux real,
                         sii_err real,
                         sii_EW_obs real,
                         siii_9069_flux real,
                         siii_9069_err real,
                         siii_9069_EW_obs real,
                         siii_9532_flux real,
                         siii_9532_err real,
                         siii_9532_EW_obs real,
                         he1_10830_flux real,
                         he1_10830_err real,
                         he1_10830_EW_obs real,
                         ContamFlag int,
                         EntryTime text
                         )''')
        self.dbConnection.commit()
        print('Done.' if self.dbCursor.rowcount > 0 else 'Table was already present.')

    def createCommentTable(self):
        print('Creating comments table in SQLLite database...')
        self.dbCursor.execute('''CREATE TABLE IF NOT EXISTS comments (
                         ParNo int,
                         ObjID int,
                         Comment text
                         )''')
        self.dbConnection.commit()
        print('Done.' if self.dbCursor.rowcount > 0 else 'Table was already present.')

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
        nonFitResults = [catalogueEntryData[-1][key] for key in catalogueEntryData[-1].keys() if key not in WISPLFDatabaseManager.fitResultKeys ]
        fitResults = { key : catalogueEntryData[-1][key] for key in WISPLFDatabaseManager.fitResultKeys }
        return tuple(nonFitResults), fitResults

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

    def saveComment(self, commentData):
        query = 'INSERT INTO comments VALUES ({})'.format(
            ','.join(['?'] * len(commentData))
        )
        self.dbCursor.execute(query, commentData)
        self.dbConnection.commit()

    def writeCatalogueTextFile(self):
        pass
