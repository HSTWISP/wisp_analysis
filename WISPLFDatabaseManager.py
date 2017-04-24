import numpy as np
import sqlite3


class WISPLFDatabaseManager:
    """
    Utility class to manage persistence of line identification data during the
    WISP line-finding process.
    """

    def __init__(self, dbFileNamePrefix):
        self.dbFileNamePrefix = dbFileNamePrefix
        self.dbFileName = '{}_sqlite.db'.format(self.dbFileNamePrefix)
        self.dbConnection = sqlite3.connect(self.dbFileName)
        self.dbCursor = self.dbConnection.cursor()
        self.checkAndInitTables()

    def __del__(self):
        self.dbCursor.commit()
        self.dbConnection.close()

    def checkAndInitTables(self):
        self.dbCursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='catalogue';")
        if self.dbCursor.fetchone() is None:
            # table does not yet exist create it
            self.createCatalogueTable()

        self.dbCursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='comments';")
        if self.dbCursor.fetchone() is None:
            # table does not yet exist create it
            self.createCommentTable()

    def createCatalogueTable(self):
        print('Creating catalogue table in SQLLite database...')
        self.dbCursor.execute('''CREATE TABLE catalogue (
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
                         ContamFlag int
                         )''')
        self.dbConnection.commit()
        print('Done.')

    def createCommentTable(self):
        print('Creating comments table in SQLLite database...')
        self.dbCursor.execute('''CREATE TABLE comments (
                         ParNo int,
                         ObjID int,
                         Comment text
                         )''')
        self.dbConnection.commit()
        print('Done.')

    def saveCatalogueEntry(self, catalogueEntryData):
        pass

    def loadCatalogueEntry(self, parNumber, objectId):
        pass

    def saveComment(self, commentData):
        pass

    def writeCatalogueTextFile(self):
        pass
