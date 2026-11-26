import unittest
from scan_world_railstatus_disturbance_guild_pairs import pack_pair,bucket
from analyze_world_railstatus_disturbance_guild_pairs import estimate,fields
class TestGuild(unittest.TestCase):
    def test_empty(self):
        self.assertIsNone(pack_pair([bucket(),bucket()],[bucket(),bucket()],[0,0],[0,0]))
        self.assertNotEqual(estimate([])["status"],"estimated")
    def test_unique_fields(self):
        for e in (False,True):
            for f in (False,True):
                q=fields(e,f);self.assertEqual(len(q),len(set(q)))
if __name__=="__main__":unittest.main()
